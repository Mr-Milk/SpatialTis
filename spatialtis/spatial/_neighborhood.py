from collections import Counter
from itertools import combinations_with_replacement, product
from typing import Callable, Mapping, Optional, Sequence
import warnings

import numpy as np
import pandas as pd
from tqdm import tqdm

from spatialtis.config import CONFIG
from spatialtis.utils import df2adata_uns, timer

from ._neighbors import Neighbors
from ._util import check_neighbors


class CellComb:
    def __init__(self, types, order):
        self.types = types
        if order:
            self.comb = [k for k in product(types, repeat=2)]
        else:
            self.comb = [k for k in combinations_with_replacement(types, 2)]
        self.relationships = {t1: [(t1, t2) for t2 in types] for t1 in types}

    def get_comb(self, t):
        return self.relationships[t]


def _count_neighbors(
        types: Sequence, relationships: Mapping, cellcomb: CellComb, storage_object: Mapping
):
    if len(relationships) > 0:
        for k, v in relationships.items():
            if len(v) > 0:
                counts = Counter([types[i] for i in v])
                center_type = types[k]
                # if center_type in counts.keys():
                for (t1, t2) in cellcomb.get_comb(center_type):
                    try:
                        c = counts[t2]
                    except (KeyError, ValueError):
                        c = 0
                    try:
                        storage_object[(t1, t2)].append(c)
                    except (KeyError, ValueError):
                        storage_object[(t2, t1)].append(c)
    itr = storage_object.keys()
    counts = [np.mean(v) if len(v) > 0 else 0 for v in storage_object.values()]
    return dict(zip(itr, counts))


def _bootstrap(
        cell_types: list, cell_neighbors: Mapping, cellcomb: CellComb, resample: int,
):
    tmp_storage = {k: [] for k in cellcomb.comb}
    real_count = _count_neighbors(cell_types, cell_neighbors, cellcomb, tmp_storage)

    perm_count = {k: [] for k in cellcomb.comb}
    shuffle_types = cell_types.copy()
    for attempt in range(resample):
        np.random.shuffle(shuffle_types)
        tmp_storage = {k: [] for k in cellcomb.comb}
        perm = _count_neighbors(shuffle_types, cell_neighbors, cellcomb, tmp_storage)
        for itr, m in perm.items():
            perm_count[itr].append(m)
    columns = cellcomb.comb
    perm_count = pd.DataFrame(perm_count, columns=columns)
    real_count = pd.DataFrame(real_count, index=[0], columns=columns)
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


def _patch_pval(
        perm_count: pd.DataFrame, real_count: pd.DataFrame, resample: int, pval: float
):
    p_gt = (perm_count >= real_count.to_numpy()).sum() / (resample + 1)
    p_lt = (perm_count <= real_count.to_numpy()).sum() / (resample + 1)
    direction = p_gt < p_lt
    redirection = [not i for i in direction]
    p = p_gt * direction + p_lt * redirection
    sig = p < pval
    sigv = (sig * np.sign(direction - 0.5)).astype(int)

    return pd.Series(sigv, index=sig.index)


def _patch_zscore(
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
        order: bool = True,
        mp: Optional[bool] = None,
):
    if mp is None:
        mp = CONFIG.MULTI_PROCESSING

    check_neighbors(n)
    cellcomb = CellComb(n.unitypes, order)
    results = dict()

    if mp & (CONFIG.OS in ["Linux", "Darwin"]):

        def exec_iterator(obj_ids):
            while obj_ids:
                done, obj_ids = ray.wait(obj_ids)
                yield ray.get(done[0])

        counts = []
        names = []
        for name, value in n.neighbors.items():
            id1 = _bootstrap_mp.remote(n.types[name], value, cellcomb, resample)
            counts.append(id1)
            names.append(name)

        for _ in tqdm(
                exec_iterator(counts),
                **CONFIG.tqdm(total=len(counts), desc="neighborhood analysis"),
        ):
            pass

        counts = ray.get(counts)

        for count, name in zip(counts, names):
            results[name] = patch_func(count[0], count[1], resample, pval)

    else:
        for name, value in tqdm(
                n.neighbors.items(), **CONFIG.tqdm(desc="neighborhood analysis"),
        ):
            [perm_count, real_count] = _bootstrap(
                n.types[name], value, cellcomb, resample
            )
            results[name] = patch_func(perm_count, real_count, resample, pval)

    results_df = pd.DataFrame(results)
    results_df.index = pd.MultiIndex.from_tuples(
        results_df.index, names=("Cell type1", "Cell type2")
    )
    results_df.rename_axis(columns=n.expobs, inplace=True)

    return results_df


def _na_fast(
        n: Neighbors,
        resample: int = 50,
        pval: float = 0.01,
        order: bool = True,
        method: str = "pval"
):
    try:
        import neighborhood_analysis as na
    except ImportError:
        return None

    check_neighbors(n)
    types = n.unitypes
    cb = na.CellCombs(types, order)

    results = {}
    for name, value in tqdm(
            n.neighbors.items(), **CONFIG.tqdm(desc="neighborhood analysis"),
    ):
        result = cb.bootstrap(n.types[name], value, resample, pval, method)
        result = {tuple(k): v for (k, v) in result}
        results[name] = result

    results_df = pd.DataFrame(results)
    results_df.index = pd.MultiIndex.from_tuples(
        results_df.index, names=("Cell type1", "Cell type2")
    )
    results_df.rename_axis(columns=n.expobs, inplace=True)

    return results_df


#@timer(prefix="Running neighborhood analysis")
def neighborhood_analysis(
        n: Neighbors,
        method: str = "pval",
        resample: int = 50,
        pval: float = 0.01,
        order: bool = True,
        export: bool = True,
        export_key: Optional[str] = None,
        return_df: bool = False,
        mp: bool = False,
):
    """Python implementation of neighborhood analysis

    Neighborhood analysis tells you the relationship between different type of cells

    There are two type of relationship, association (1) or avoidance (-1), no relationship (0).

    Args:
        n: A spatialtis.Neighbors object, neighbors are already computed
        method: "pval" or "zscore"
        resample: perform resample for how many times
        pval: if smaller than pval, reject null hypothesis (No relationship)
        order: if False, Cell A - Cell B and Cell B - Cell A are the same interaction.
        export: whether to export the result to anndata.uns
        export_key: the key used to export
        return_df: whether to return the result
        mp: whether to enable multiprocessing

    .. seealso:: `spatial_enrichment_analysis <#spatialtis.plotting.spatial_enrichment_analysis>`_


    """

    if export_key is None:
        export_key = CONFIG.neighborhood_analysis_key
    else:
        CONFIG.neighborhood_analysis_key = export_key

    try:
        df = _na_fast(n, resample, pval, order, method)
    except ImportError:
        warnings.warn("Try `pip install neighborhood_analysis` which is much faster than current one.")
        if method == "pval":
            df = _main(n, _patch_pval, resample=resample, pval=pval, order=order, mp=mp)

        else:
            df = _main(n, _patch_zscore, resample=resample, pval=pval, order=order, mp=mp)

    if method == "pval":
        df = df.T.astype(int)
    else:
        df = df.T

    if export:
        df2adata_uns(df, n.adata, export_key, params={"order": order, "method": method, "pval": pval})

    if return_df:
        return df


@timer(prefix="Running spatial enrichment analysis")
def spatial_enrichment_analysis(
        n: Neighbors,
        threshold: Optional[float] = None,
        layers_key: Optional[str] = None,
        resample: int = 50,
        export: bool = True,
        export_key: Optional[str] = None,
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
            export: whether to export the result to anndata.uns
            export_key: the key used to export
            return_df: whether to return the result
            mp: whether to enable multiprocessing

        .. seealso:: `neighborhood_analysis <#spatialtis.plotting.neighborhood_analysis>`_

        """

    check_neighbors(n)
    data = n.adata

    if export_key is None:
        export_key = CONFIG.spatial_enrichment_analysis_key
    else:
        CONFIG.spatial_enrichment_analysis_key = export_key

    if (threshold is None) & (layers_key is None):
        raise ValueError("Either specific a threshold or a layers key.")
    elif layers_key is not None:
        if threshold is not None:
            warnings.warn("You specific both threshold and layers_key, using user defined layers_key")
        CONFIG.spatial_enrichment_analysis_layers_key = layers_key
    elif threshold is not None:
        layers_key = CONFIG.spatial_enrichment_analysis_layers_key
        data[layers_key] = (data >= threshold)





    df = _main(n, _patch_zscore, resample=resample, mp=mp)

    df = df.T

    if export:
        df2adata_uns(df, n.adata, export_key)

    if return_df:
        return df
