from collections import Counter
from typing import Optional, Sequence, Union

import numpy as np
import pandas as pd
from anndata import AnnData

from ..utils import df2adata_uns

Num = Union[int, float]


def type_counter(
    adata: AnnData,
    groupby: Union[Sequence, str],
    type_col: str,
    selected_types: Optional[Sequence] = None,
) -> pd.DataFrame:

    if isinstance(groupby, str):
        groupby = [groupby]
    else:
        groupby = list(groupby)
    df = adata.obs[groupby + [type_col]]

    if selected_types is not None:
        df = df[df[type_col].isin(selected_types)]

    groups = df.groupby(groupby)
    types = pd.unique(df[type_col])
    # TODO: sorted types based on selected types
    matrix = list()
    mindex = list()
    for i, (n, g) in enumerate(groups):
        c = Counter(g[type_col])
        matrix.append([c[t] for t in types])
        # one column index
        if isinstance(n, str):
            mindex.append((n, i,))
        # multiIndex
        else:
            mindex.append((*n, i,))

    multiIndex = pd.MultiIndex.from_tuples(mindex)
    all_results = pd.DataFrame(dict(zip(types, np.array(matrix).T)), index=multiIndex)
    all_results.index.set_names(groupby + ["id"], inplace=True)
    all_results.columns.set_names("type", inplace=True)

    return all_results  # obs as index, var as columns


def cell_components(
    adata: AnnData,
    type_col: str,
    groupby: list,
    selected_types: Optional[Sequence] = None,
    export: bool = True,
    export_key: str = "cell_components",
    return_df: bool = False,
):
    results = type_counter(adata, groupby, type_col, selected_types=selected_types)

    if export:
        df2adata_uns(results, adata, export_key)

    if return_df:
        return results


def cell_co_occurrence(
    adata: AnnData,
    type_col: str,
    groupby: list,
    selected_types: list = None,
    export: bool = True,
    export_key: str = "co_occurrence",
    threshold: int = 50,
    return_df: bool = False,
):
    counter = type_counter(adata, groupby, type_col, selected_types=selected_types)

    occurrence = counter.gt(threshold)

    if export:
        df2adata_uns(occurrence, adata, export_key)

    if return_df:
        return occurrence


def cell_density(
    adata: AnnData,
    type_col: str,
    groupby: list,
    size: Union[Sequence[Sequence[Num]], Sequence[Num]],
    selected_types: Optional[Sequence] = None,
    export: bool = True,
    export_key: str = "cell_density",
    ratio: float = 1.0,
    return_df: bool = False,
):
    counter = type_counter(adata, groupby, type_col, selected_types=selected_types)
    if isinstance(size[0], (int, float)):
        area = size[0] * size[1]
        results = counter.div(area * ratio)

    elif isinstance(size[0], Sequence):
        area = pd.Series({"area": [s[0] * s[1] * ratio for s in size]})
        results = counter.div(area["area"], axis=0)
    else:
        raise ValueError("Unrecognized size input")

    if export:
        df2adata_uns(results, adata, export_key)

    if return_df:
        return results

def cell_morphology(
    adata: AnnData,
    type_col: str,
    groupby: list,
    selected_types: Optional[Sequence] = None,
    area_col: str = "area",
    eccentricity_col: str = "eccentricity",
    export: bool = True,
    export_key: str = "cell_components",
    return_df: bool = False,
):
    pass