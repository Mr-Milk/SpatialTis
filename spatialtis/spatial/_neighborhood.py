from collections import Counter
from itertools import product
from typing import Callable, Mapping, Sequence

import numpy as np
import pandas as pd

from spatialtis.utils import df2adata_uns
from ._neighbors import Neighbors
from ._util import check_neighbors


def _count_neighbors(
        types: Sequence,
        relationships: Mapping,
        storage_object: Mapping
):
    if len(relationships) > 0:
        for k, v in relationships.items():
            if len(v) > 0:
                for t, count in Counter([types[i] for i in v]).items():
                    select = (types[k], t,)
                    storage_object[select].append(count)

    itr = storage_object.keys()
    counts = [np.mean(v) if len(v) > 0 else 0 for v in storage_object.values()]
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
    for attempt in range(0, resample):
        np.random.shuffle(shuffle_types)
        tmp_storage = {k: [] for k in cell_interactions}
        perm = _count_neighbors(shuffle_types, cell_neighbors, tmp_storage)
        for itr, m in perm.items():
            perm_count[itr].append(m)
    perm_count = pd.DataFrame(perm_count, columns=cell_interactions)
    real_count = pd.DataFrame(real_count, index=[0], columns=cell_interactions)

    return perm_count, real_count


def _patch_neighborhood(
        perm_count: pd.DataFrame,
        real_count: pd.DataFrame,
        resample: int,
        pval: float
):
    p_gt = (perm_count >= real_count.to_numpy()).sum() / (resample + 1)
    p_lt = (perm_count <= real_count.to_numpy()).sum() / (resample + 1)
    direction = p_lt < p_gt
    p = p_lt * direction + p_gt * (direction is False)
    sig = p < pval
    sigv = sig * np.sign((0.5 - direction))
    return sigv


def _patch_spatial_enrichment(
        perm_count: pd.DataFrame,
        real_count: pd.DataFrame,
        *args
):

    z_score = (real_count.to_numpy()[0] - perm_count.mean()) / perm_count.std()
    # inf ---> NaN, and then NaN ---> 0
    z_score = z_score.replace([np.inf, -np.inf], np.nan).fillna(0)

    return z_score


def _main(
        n: Neighbors,
        patch_func: Callable,
        resample: int = 50,
        pval: float = 0.01
):
    check_neighbors(n)
    cell_interactions = [k for k in product(n.unitypes, repeat=2)]
    results = dict()

    for name, value in n.neighbors.items():
        perm_count, real_count = _bootstrap(n.types[name], value, cell_interactions, resample)
        results[name] = patch_func(perm_count, real_count, resample, pval)

    results_df = pd.DataFrame(results)
    results_df.index = pd.MultiIndex.from_tuples(results_df.index, names=('Cell type1', 'Cell type2'))
    results_df.rename_axis(columns=n.expobs, inplace=True)

    return results_df


def neighborhood_analysis(
        n: Neighbors,
        resample: int = 50,
        pval: float = 0.01,
        export: bool = True,
        export_key: str = "neighborhood_analysis",
        return_df: bool = False,
        overwrite: bool = False,
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
        overwrite: whether to overwrite your previous results (if existed)

    .. sea also:: :function: `spatial_enrichment_analysis`


    """
    df = _main(n, _patch_neighborhood, resample=resample, pval=pval)

    df = df.T.astype(int)

    if export:
        df2adata_uns(df, n.adata, export_key, overwrite)

    if return_df:
        return df


def spatial_enrichment_analysis(
        n: Neighbors,
        resample: int = 50,
        export: bool = True,
        export_key: str = "spatial_enrichment_analysis",
        return_df: bool = False,
        overwrite: bool = False,
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
            overwrite: whether to overwrite your previous results (if existed)

        .. sea also:: :function: `neighborhood_analysis`

        """
    df = _main(n, _patch_spatial_enrichment, resample=resample)

    df = df.T

    if export:
        df2adata_uns(df, n.adata, export_key, overwrite)

    if return_df:
        return df
