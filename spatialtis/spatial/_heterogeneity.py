import warnings
from typing import Optional, Sequence, Union

import pandas as pd
from anndata import AnnData
from scipy.stats import entropy
from spatialentropy import altieri_entropy, leibovici_entropy
from tqdm import tqdm

from spatialtis.config import CONFIG

from ..sta.statistics import type_counter
from ..utils import df2adata_uns, lprint, timer

if CONFIG.OS in ["Linux", "Darwin"]:
    try:
        import ray
    except ImportError:
        raise ImportError(
            "You don't have ray installed or your OS don't support ray.",
            "Try `pip install ray` or use `mp=False`",
        )

    @ray.remote
    def altieri_entropy_mp(*args, **kwargs):
        e = altieri_entropy(*args, **kwargs)
        return e.entropy

    @ray.remote
    def leibovici_entropy_mp(*args, **kwargs):
        e = leibovici_entropy(*args, **kwargs)
        return e.entropy


@timer(prefix="Running plotting heterogeneity")
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
) -> Optional[pd.DataFrame]:
    """compute plotting heterogeneity

    Here we use entropy for plotting heterogeneity, which describes the amount of information.
    To compare the difference within a group (eg. different samples from same tumor), Kullback–Leibler divergences
    for each sample within the group are computed, smaller value indicates less difference within group.

    Args:
        adata: anndata object to perform analysis
        groupby: list of names describes your experiments
        method: Use which entropy: "shannon", "leibovici", "altieri"
        base: the log base
        d: the distance threshold to determine co-occurrence events (method="leibovici")
        cut: distance interval (method="altieri")
        compare: Compute Kullback-Leibler divergences based on which level
        type_key: the key of cell type in anndata.obs (Default: spatialtis.CONFIG.CELL_TYPE_KEY)
        centroid_key: the key of cell centroid in anndata.obs (Default: spatialtis.CONFIG.CENTROID_KEY)
        export: whether to export the result to anndata.uns
        export_key: the key used to export
        return_df: whether to return the result
        mp: whether to enable multiprocessing

    """
    if groupby is None:
        groupby = CONFIG.EXP_OBS
    if type_key is None:
        type_key = CONFIG.CELL_TYPE_KEY
    if centroid_key is None:
        centroid_key = CONFIG.CENTROID_KEY

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
        lprint("Method: Shannon entropy")
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
            lprint("Method: Altieri entropy")
        if method == "leibovici":
            lprint("Method: Leibovici entropy")

        df = adata.obs[groupby + [type_key, centroid_key]]

        ent = list()
        mindex = list()
        groups = df.groupby(groupby)

        if mp & (CONFIG.OS in ["Linux", "Darwin"]):

            def exec_iterator(obj_ids):
                while obj_ids:
                    done, obj_ids = ray.wait(obj_ids)
                    yield ray.get(done[0])

            results = list()

            for i, (n, g) in enumerate(groups):
                types = list(g[type_key])
                points = [eval(i) for i in g[centroid_key]]
                if method == "altieri":
                    results.append(
                        altieri_entropy_mp.remote(points, types, cut=cut, base=base)
                    )
                else:
                    results.append(
                        leibovici_entropy_mp.remote(points, types, d=d, base=base)
                    )
                # one column index
                if isinstance(n, str):
                    mindex.append((n, i,))
                # multiIndex
                else:
                    mindex.append((*n, i,))

            for _ in tqdm(
                exec_iterator(results),
                **CONFIG.tqdm(total=len(results), desc="heterogeneity")
            ):
                pass

            mp_results = ray.get(results)

            for e in mp_results:
                ent.append(e)

        else:
            for i, (n, g) in enumerate(
                tqdm(groups, **CONFIG.tqdm(desc="heterogeneity",))
            ):
                types = list(g[type_key])
                points = [eval(i) for i in g[centroid_key]]
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
