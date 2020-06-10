from typing import Optional, Sequence, Union

import numpy as np
import pandas as pd
from anndata import AnnData
from scipy.stats import entropy
from spatialentropy import altieri_entropy, leibovici_entropy
from tqdm import tqdm

from spatialtis.config import CONFIG

from ..sta.statistics import type_counter
from ..utils import df2adata_uns

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


def spatial_heterogeneity(
    adata: AnnData,
    groupby: Union[Sequence, str, None] = None,
    type_col: Optional[str] = None,
    centroid_col: Optional[str] = None,
    method: str = "altieri",  # shannon, leibovici, altieri
    base: Union[int, float, None] = None,
    d: Optional[int] = None,
    cut: Union[int, Sequence, None] = None,
    compare: Optional[str] = None,
    export_key: str = "spatial_heterogeneity",
    return_df: bool = False,
    mp: Optional[bool] = None,
) -> Optional[pd.DataFrame]:
    """compute spatial heterogeneity
    Here we use entropy for spatial heterogeneity, which describes the amount of information.
    To compare the difference within a group (eg. different samples from same tumor), Kullback–Leibler divergences
    for each sample within the group are computed, smaller value indicates less difference within group.

    Args:
        adata: anndata object to perform analysis
        groupby: list of names describes your experiments
        type_col: name of the cell types column
        centroid_col:
        method:
        base:
        d:
        cut:
        compare: Compute Kullback-Leibler divergences based on which level
        export_key: the key name to store info, exported to anndata.uns field
        return_df: whether to return a pandas DataFrame object
        mp:

    """
    if groupby is None:
        groupby = CONFIG.EXP_OBS
    if type_col is None:
        type_col = CONFIG.CELL_TYPE_COL
    if centroid_col is None:
        centroid_col = CONFIG.CENTROID_COL
    if mp is None:
        mp = CONFIG.MULTI_PROCESSING

    if method not in ["shannon", "altieri", "leibovici"]:
        raise ValueError(
            "Available entropy methods are 'shannon', 'altieri', 'leibovici'."
        )

    if method == "shannon":
        df = type_counter(adata, groupby, type_col)

        if len(df.columns) == 1:
            print("No heterogeneity, only one type of cell found.")
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
        df = adata.obs[groupby + [type_col, centroid_col]]

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
                types = list(g[type_col])
                points = [eval(i) for i in g[centroid_col]]
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
                total=len(results),
                desc="heterogeneity",
                bar_format=CONFIG.PBAR_FORMAT,
                disable=(not CONFIG.PROGRESS_BAR),
            ):
                pass

            mp_results = ray.get(results)

            for e in mp_results:
                ent.append(e)

        else:
            for i, (n, g) in enumerate(
                tqdm(
                    groups,
                    desc="heterogeneity",
                    bar_format=CONFIG.PBAR_FORMAT,
                    disable=(not CONFIG.PROGRESS_BAR),
                )
            ):
                types = list(g[type_col])
                points = [eval(i) for i in g[centroid_col]]
                if method == "altieri":
                    e = altieri_entropy(points, types, cut=cut, base=base)
                else:
                    e = leibovici_entropy(points, types, d=d, base=base)
                ent.append(e.entropy)
                mindex.append(n)
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
    df2adata_uns(roi_heterogeneity, adata, export_key)

    if return_df:
        return roi_heterogeneity
