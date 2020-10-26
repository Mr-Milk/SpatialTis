import warnings
from ast import literal_eval
from typing import Optional, Sequence, Union

import pandas as pd
from anndata import AnnData
from scipy.stats import entropy
from spatialentropy import altieri_entropy, leibovici_entropy
from tqdm import tqdm

from spatialtis.config import CONFIG
from spatialtis.stats.statistics import type_counter
from spatialtis.utils import (
    create_remote,
    df2adata_uns,
    get_default_params,
    log_print,
    reuse_docstring,
    run_ray,
    timer,
)


@timer(prefix="Running spatial heterogeneity")
@get_default_params
@reuse_docstring()
def spatial_heterogeneity(
    adata: AnnData,
    groupby: Union[Sequence, str, None] = None,
    method: str = "leibovici",  # shannon, leibovici, altieri
    base: Union[int, float, None] = None,
    d: Optional[int] = None,
    cut: Union[int, Sequence, None] = None,
    compare: Optional[str] = None,
    type_key: Optional[str] = None,
    centroid_key: Optional[str] = None,
    export: bool = True,
    export_key: Optional[str] = None,
    return_df: bool = False,
    mp: Optional[bool] = None,
):
    """Evaluate tissue heterogeneity based on entropy

    Entropy describes the amount of information.

    - `Shannon entropy <../about/implementation.html#shannon-entropy>`_ (No spatial info included):\
        To compare the difference within a group (eg. different samples from same tumor), Kullback–Leibler divergences\
        for each sample within the group are computed, smaller value indicates less difference within group.
    - `Leibovici entropy <../about/implementation.html#leibovici-entropy>`_:\
    You can specific the distance threshold to determine co-occurrence events.
    - `Altieri entropy <../about/implementation.html#altieri-entropy>`_:\
    You can specific the distance interval to determine co-occurrence events.

    Args:
        adata: {adata}
        groupby: {groupby}
        method: Options are "shannon", "leibovici" and "altieri"
        base: The log base
        d: The distance threshold to determine co-occurrence events (method="leibovici")
        cut: Distance interval (method="altieri")
        compare: Compute Kullback-Leibler divergences based on which level
        type_key: {type_key}
        centroid_key: {centroid_key}
        export: {export}
        export_key: {export_key}
        return_df: {return_df}
        mp: {mp}

    """
    if export_key is None:
        export_key = CONFIG.spatial_heterogeneity_key
    else:
        CONFIG.spatial_heterogeneity_key = export_key

    if mp is None:
        mp = CONFIG.MULTI_PROCESSING

    if method not in ["shannon", "altieri", "leibovici"]:
        raise ValueError(
            "Available entropy methods are 'shannon', 'altieri', 'leibovici'."
        )

    if method == "shannon":
        log_print("Method: Shannon entropy")
        df = type_counter(adata, groupby, type_key=type_key)

        if len(df.columns) == 1:
            warnings.warn(
                "No heterogeneity, you only have one type of cell.", UserWarning
            )
            return None

        KL_div = dict()
        if compare is not None:
            compare_index = CONFIG.EXP_OBS.index(compare)
            groups = df.groupby(level=compare)
            for n, g in groups:
                count_g = g.sum()
                qk = count_g.div(count_g.sum())
                KL_div[n] = qk

        ent, KL, KL_level = list(), list(), list()
        for row in df.iterrows():
            pk = list(row[1].div(row[1].sum()))
            if compare is not None:
                compare_level = row[0][compare_index]
                KL.append(entropy(pk, KL_div[compare_level], base=base))
                KL_level.append(compare_level)
            ent.append(entropy(pk, base=base))

        data = {"heterogeneity": ent}
        if compare is not None:
            data["KL"] = KL
            data["level"] = KL_level
        roi_heterogeneity = pd.DataFrame(data=data, index=df.index)

    else:
        if method == "altieri":
            log_print("Method: Altieri entropy")
        if method == "leibovici":
            log_print("Method: Leibovici entropy")

        df = adata.obs[groupby + [type_key, centroid_key]]

        ent = list()
        mindex = list()
        groups = df.groupby(groupby)

        if mp:

            def altieri_entropy_mp(*args, **kwargs):
                return altieri_entropy(*args, **kwargs).entropy

            def leibovici_entropy_mp(*args, **kwargs):
                return leibovici_entropy(*args, **kwargs).entropy

            altieri_entropy_mp, leibovici_entropy_mp = create_remote(
                [altieri_entropy_mp, leibovici_entropy_mp]
            )

            jobs = []
            for i, (n, g) in enumerate(groups):
                types = list(g[type_key])
                points = [literal_eval(i) for i in g[centroid_key]]
                if method == "altieri":
                    jobs.append(
                        altieri_entropy_mp.remote(points, types, cut=cut, base=base)
                    )
                else:
                    jobs.append(
                        leibovici_entropy_mp.remote(points, types, d=d, base=base)
                    )
                # one column index
                if isinstance(n, str):
                    mindex.append((n, i,))
                # multiIndex
                else:
                    mindex.append((*n, i,))

            mp_results = run_ray(
                jobs,
                tqdm_config=CONFIG.tqdm(
                    total=len(jobs), desc="Calculating heterogeneity"
                ),
            )

            for e in mp_results:
                ent.append(e)

        else:
            for i, (n, g) in enumerate(
                tqdm(groups, **CONFIG.tqdm(desc="Calculating heterogeneity",))
            ):
                types = list(g[type_key])
                points = [literal_eval(i) for i in g[centroid_key]]
                if method == "altieri":
                    e = altieri_entropy(points, types, cut=cut, base=base)
                else:
                    e = leibovici_entropy(points, types, d=d, base=base)
                ent.append(e.entropy)
                # one column index
                if isinstance(n, str):
                    mindex.append((n, i,))
                # multiIndex
                else:
                    mindex.append((*n, i,))
        data = {"heterogeneity": ent}
        roi_heterogeneity = pd.DataFrame(data=data)
        roi_heterogeneity.index = pd.MultiIndex.from_tuples(
            mindex, names=groupby + ["id"]
        )

    # export to anndata
    if export:
        df2adata_uns(roi_heterogeneity, adata, export_key, params={"method": method})

    if return_df:
        return roi_heterogeneity
