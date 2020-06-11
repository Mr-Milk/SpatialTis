from collections import Counter
from itertools import product
from typing import Callable, Mapping, Optional, Sequence

import numpy as np
import pandas as pd
from tqdm import tqdm

from spatialtis.config import CONFIG
from spatialtis.utils import df2adata_uns

from ._neighbors import Neighbors
from ._util import check_neighbors


def _count_neighbors(types: Sequence, relationships: Mapping, storage_object: Mapping):
    cells = len(relationships)
    if cells > 0:
        for k, v in relationships.items():
            if len(v) > 0:
                for t, count in Counter([types[i] for i in v]).items():
                    select = (types[k], t)
                    storage_object[select].append(count)

    itr = storage_object.keys()
    counts = [(np.sum(v) / cells) for v in storage_object.values()]
    return dict(zip(itr, counts))


def _bootstrap(
    cell_types: list,
    cell_neighbors: Mapping,
    cell_interactions: Sequence,
    resample: int,
):
    tmp_storage = {k: [] for k in cell_interactions}
    real_count = _count_neighbors(cell_types, cell_neighbors, tmp_storage)

    perm_count = {k: [] for k in cell_interactions}
    shuffle_types = cell_types.copy()
    for attempt in range(resample):
        np.random.shuffle(shuffle_types)
        tmp_storage = {k: [] for k in cell_interactions}
        perm = _count_neighbors(shuffle_types, cell_neighbors, tmp_storage)
        for itr, m in perm.items():
            perm_count[itr].append(m)
    perm_count = pd.DataFrame(perm_count, columns=cell_interactions)
    real_count = pd.DataFrame(real_count, index=[0], columns=cell_interactions)
    return [perm_count, real_count]


if CONFIG.OS in ["Linux", "Darwin"]:
    try:
        import ray
    except ImportError:
        raise ImportError(
            "You don't have ray installed or your OS don't support ray.",
            "Try `pip install ray` or use `mp=False`",
        )

    _bootstrap_mp = ray.remote(_bootstrap)


def _patch_neighborhood(
    perm_count: pd.DataFrame, real_count: pd.DataFrame, resample: int, pval: float
):
    p_gt = (perm_count >= real_count.to_numpy()).sum() / (resample + 1)
    p_lt = (perm_count <= real_count.to_numpy()).sum() / (resample + 1)
    direction = p_lt < p_gt
    re_direction = [not i for i in direction]
    p = p_gt * direction + p_lt * re_direction
    sig = (1 - p) < pval

    sigv = []
    for pv, s, d in zip(p, sig, direction):
        if not s:
            sigv.append(0)
        else:
            if d:
                sigv.append(1)
            else:
                sigv.append(-1)

    return pd.Series(sigv, index=sig.index)


def _patch_spatial_enrichment(
    perm_count: pd.DataFrame, real_count: pd.DataFrame, *args
):
    z_score = (real_count.to_numpy()[0] - perm_count.mean()) / perm_count.std()
    # inf ---> NaN, and then NaN ---> 0
    z_score = z_score.replace([np.inf, -np.inf], np.nan).fillna(0)

    return z_score


def _main(
    n: Neighbors,
    patch_func: Callable,
    resample: int = 50,
    pval: float = 0.01,
    mp: Optional[bool] = None,
):
    if mp is None:
        mp = CONFIG.MULTI_PROCESSING

    check_neighbors(n)
    cell_interactions = [k for k in product(n.unitypes, repeat=2)]
    results = dict()

    if mp & (CONFIG.OS in ["Linux", "Darwin"]):

        def exec_iterator(obj_ids):
            while obj_ids:
                done, obj_ids = ray.wait(obj_ids)
                yield ray.get(done[0])

        counts = []
        names = []
        for name, value in n.neighbors.items():
            id1 = _bootstrap_mp.remote(
                n.types[name], value, cell_interactions, resample
            )
            counts.append(id1)
            names.append(name)

        for _ in tqdm(
            exec_iterator(counts),
            total=len(counts),
            unit="ROI",
            desc="neighborhood analysis",
            bar_format=CONFIG.PBAR_FORMAT,
            disable=(not CONFIG.PROGRESS_BAR),
        ):
            pass

        counts = ray.get(counts)

        for count, name in zip(counts, names):
            results[name] = patch_func(count[0], count[1], resample, pval)

    else:
        for name, value in tqdm(
            n.neighbors.items(),
            unit="ROI",
            desc="neighborhood analysis",
            bar_format=CONFIG.PBAR_FORMAT,
            disable=(not CONFIG.PROGRESS_BAR),
        ):
            [perm_count, real_count] = _bootstrap(
                n.types[name], value, cell_interactions, resample
            )
            results[name] = patch_func(perm_count, real_count, resample, pval)

    results_df = pd.DataFrame(results)
    results_df.index = pd.MultiIndex.from_tuples(
        results_df.index, names=("Cell type1", "Cell type2")
    )
    results_df.rename_axis(columns=n.expobs, inplace=True)

    return results_df


def neighborhood_analysis(
    n: Neighbors,
    resample: int = 50,
    pval: float = 0.01,
    export: bool = True,
    export_key: str = "neighborhood_analysis",
    return_df: bool = False,
    mp: bool = False,
):
    """Python implementation of histocat's neighborhood analysis

    Neighborhood analysis tells you the relationship between different type of cells

    There are two type of relationship, association (1) or avoidance (-1), no relationship (0).

    Args:
        n: A spatialtis.Neighbors object, neighbors are already computed
        resample: perform resample for how many times
        pval: if smaller than pval, reject null hypothesis (No relationship)
        export: whether export to anndata object uns field
        export_key: which key used to export
        return_df: whether to return result dataframe
        mp: whether enable parallel processing

    .. seealso:: `spatial_enrichment_analysis <#spatialtis.spatial.spatial_enrichment_analysis>`_


    """
    df = _main(n, _patch_neighborhood, resample=resample, pval=pval, mp=mp)

    df = df.T.astype(int)

    if export:
        df2adata_uns(df, n.adata, export_key)

    if return_df:
        return df


def spatial_enrichment_analysis(
    n: Neighbors,
    resample: int = 50,
    export: bool = True,
    export_key: str = "spatial_enrichment_analysis",
    return_df: bool = False,
    mp: bool = False,
):
    """An alternative neighborhood analysis

        Neighborhood analysis tells you the relationship between different type of cells

        This method is purposed in MIBI's paper, the major difference is that this method used z-score for relationships

        To be noticed, when there is no cells, this will cause NaN, those value will be replaced with 0

        Args:
            n: A spatialtis.Neighbors object, neighbors are already computed
            resample: perform resample for how many times
            export: whether export to anndata object uns field
            export_key: which key used to export
            return_df: whether to return result dataframe
            mp: whether enable parallel processing

        .. seealso:: `neighborhood_analysis <#spatialtis.spatial.neighborhood_analysis>`_

        """
    df = _main(n, _patch_spatial_enrichment, resample=resample, mp=mp)

    df = df.T

    if export:
        df2adata_uns(df, n.adata, export_key)

    if return_df:
        return df
