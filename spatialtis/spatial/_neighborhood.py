import numpy as np
import pandas as pd

from collections import Counter
from itertools import product
from typing import Sequence, Mapping, Callable

from ._util import check_neighbors

from ._neighbors import Neighbors


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
        cell_types: Sequence,
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
        perm = _count_neighbors(cell_types, cell_neighbors, tmp_storage)
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
    z_score = (real_count.to_numpy() - perm_count.mean()) / perm_count.std()
    return z_score


def _main(
        n: Neighbors,
        patch_func: Callable,
        resample: int = 50,
        pval: float = 0.01
):
    check_neighbors(n)
    cell_interactions = [k for k in product(n.types, repeat=2)]
    results = dict()

    for name, value in n.neighbors.items():
        perm_count, real_count = _bootstrap(n.types[name], value, cell_interactions, resample)
        patch_func(perm_count, real_count, resample, pval)

    results_df = pd.DataFrame(results)
    results_df.index = pd.MultiIndex.from_tuples(results_df.index, names=('type1', 'type2'))

    return results_df


def neighborhood_analysis(
        n: Neighbors,
        resample: int = 50,
        pval: float = 0.01
):
    """Python implementation of histocat's neighborhood analysis

    Neighborhood analysis tells you the relationship between different type of cells
    There are two type of relationship, association (1) or avoidance (-1), no relationship (0).

    Args:
        n: A spatialtis.Neighbors object, neighbors are already computed
        resample: perform resample for how many times
        pval: if smaller than pval, reject null hypothesis (No relationship)

    .. sea also:: :function: `spatial_enrichment_analysis`


    """
    df = _main(n, _patch_neighborhood, resample=resample, pval=pval)

    return df.T.astype(int)


def spatial_enrichment_analysis(
        n: Neighbors,
        resample: int = 50
):
    """An alternative neighborhood analysis

        Neighborhood analysis tells you the relationship between different type of cells
        This method is purposed in MIBI's paper, the major difference is that this method used z-score for relationships

        Args:
            n: A spatialtis.Neighbors object, neighbors are already computed
            resample: perform resample for how many times

        .. sea also:: :function: `neighborhood_analysis`

        """
    df = _main(n, _patch_spatial_enrichment, resample=resample)

    return df.T
