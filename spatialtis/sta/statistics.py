from collections import Counter
from typing import Optional, Sequence, Union

import numpy as np
import pandas as pd
from anndata import AnnData

from spatialtis.config import CONFIG
from spatialtis.utils import df2adata_uns

Num = Union[int, float]


def type_counter(
    adata: AnnData,
    groupby: Union[Sequence, str, None] = None,
    type_col: Optional[str] = None,
    selected_types: Optional[Sequence] = None,
) -> pd.DataFrame:
    """(private) To count how many type of value in a col

    Args:
        adata: anndata object to perform analysis
        groupby: how your experiments grouped, (Default: read from spatialtis.CONFIG.EXP_OBS)
        type_col: the key name of cell type in anndata.obs (Default: read from spatialtis.CONFIG.CELL_TYPE_COL)
        selected_types: selected cell types you want to count

    """
    if groupby is None:
        if CONFIG.EXP_OBS is not None:
            groupby = CONFIG.EXP_OBS
        else:
            raise ValueError(
                "Experiment observation unclear, set `spatialtis.CONFIG.EXP_OBS` or use argument "
                "`groupby=`"
            )
    else:
        groupby = list(groupby)

    if type_col is None:
        type_col = CONFIG.CELL_TYPE_COL

    df = adata.obs[groupby + [type_col]]

    if selected_types is not None:
        df = df[df[type_col].isin(selected_types)]
    # the order of type will follow the order in df
    groups = df.groupby(groupby)
    types = pd.unique(df[type_col])
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
    groupby: Union[Sequence, str, None] = None,
    type_col: Optional[str] = None,
    selected_types: Optional[Sequence] = None,
    export: bool = True,
    export_key: str = "cell_components",
    return_df: bool = False,
    overwrite: bool = False,
):
    """the proportion of different type of cells

    Args:
        adata: anndata object to perform analysis
        groupby: how your experiments grouped, (Default: read from spatialtis.CONFIG.EXP_OBS)
        type_col: the key name of cell type in anndata.obs (Default: read from spatialtis.CONFIG.CELL_TYPE_COL)
        selected_types: selected cell types you want to count
        export: whether to export to anndata.uns field
        export_key: the key name that used to record the results in anndata.uns field (Default: "cell_components")
        return_df: whether to return an pandas.DataFrame
        overwrite: whether to overwrite if the key existed

    Return:
        pandas.DataFrame

    """

    counter = type_counter(adata, groupby, type_col, selected_types)

    if export:
        df2adata_uns(counter, adata, export_key, overwrite)

    if return_df:
        return counter


def cell_co_occurrence(
    adata: AnnData,
    groupby: Union[Sequence, str, None] = None,
    type_col: Optional[str] = None,
    selected_types: Optional[Sequence] = None,
    export: bool = True,
    export_key: str = "cell_co_occurrence",
    threshold: int = 50,
    return_df: bool = False,
    overwrite: bool = False,
):
    """The probability of two type of cells occur simultaneously

        Args:
            adata: anndata object to perform analysis
            groupby: how your experiments grouped, (Default: read from spatialtis.CONFIG.EXP_OBS)
            type_col: the key name of cell type in anndata.obs (Default: read from spatialtis.CONFIG.CELL_TYPE_COL)
            selected_types: selected cell types you want to count
            export: whether to export to anndata.uns field
            export_key: the key name that used to record the results in anndata.uns field (Default: "cell_co_occurrence")
            threshold: this value determines the presence/absence of a cell type in ROI
            return_df: whether to return an pandas.DataFrame
            overwrite: whether to overwrite if the key existed

        Return:
            pandas.DataFrame

        """
    counter = type_counter(adata, groupby, type_col, selected_types)

    occurrence = counter.gt(threshold)

    if export:
        df2adata_uns(occurrence, adata, export_key, overwrite)

    if return_df:
        return occurrence


def cell_density(
    adata: AnnData,
    size: Union[Sequence[Sequence[Num]], Sequence[Num]],
    ratio: float = 1.0,
    groupby: Union[Sequence, str, None] = None,
    type_col: Optional[str] = None,
    selected_types: Optional[Sequence] = None,
    export: bool = True,
    export_key: str = "cell_density",
    return_df: bool = False,
    overwrite: bool = False,
):
    """Calculating cell density in each ROI

    Args:
        adata: anndata object to perform analysis
        size: the size of each ROI, if some ROI is different, please specific each ROI size in a list, the order must
        follow the index order in anndata object.
        ratio: the ratio between pixel and real size, default is 1.0
        groupby: how your experiments grouped, (Default: read from spatialtis.CONFIG.EXP_OBS)
        type_col: the key name of cell type in anndata.obs (Default: read from spatialtis.CONFIG.CELL_TYPE_COL)
        selected_types: selected cell types you want to count
        export: whether to export to anndata.uns field
        export_key: the key name that used to record the results in anndata.uns field (Default: "cell_density")
        return_df: whether to return an pandas.DataFrame
        overwrite: whether to overwrite if the key existed

    Return:
            pandas.DataFrame

    """
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
        df2adata_uns(results, adata, export_key, overwrite)

    if return_df:
        return results


def cell_morphology(
    adata: AnnData,
    metrics_col: str = "eccentricity",
    groupby: Union[Sequence, str, None] = None,
    type_col: Optional[str] = None,
    selected_types: Optional[Sequence] = None,
    export: bool = True,
    export_key: str = "cell_morphology",
    return_df: bool = False,
    overwrite: bool = False,
):
    """Cell morphology variation between different groups

        Args:
            adata: anndata object to perform analysis
            metrics_col: which col was used to measure cell morphology
            (Default: "eccentricity", Options: "eccentricity", "area")
            groupby: how your experiments grouped, (Default: read from spatialtis.CONFIG.EXP_OBS)
            type_col: the key name of cell type in anndata.obs (Default: read from spatialtis.CONFIG.CELL_TYPE_COL)
            selected_types: selected cell types you want to count
            export: whether to export to anndata.uns field
            export_key: the key name that used to record the results in anndata.uns field (Default: "cell_morphology")
            return_df: whether to return an pandas.DataFrame
            overwrite: whether to overwrite if the key existed

        Return:
                pandas.DataFrame

    """

    if groupby is None:
        groupby = CONFIG.EXP_OBS
    else:
        groupby = list(groupby)

    if type_col is None:
        type_col = CONFIG.CELL_TYPE_COL

    key = groupby + [type_col, metrics_col]
    df = adata.obs[key]

    if selected_types is not None:
        df = df[df[type_col].isin(selected_types)]

    df.index = pd.MultiIndex.from_frame(df[groupby + [type_col]])
    df = df.drop(groupby + [type_col], axis=1)

    if export:
        df2adata_uns(df, adata, export_key, overwrite)

    if return_df:
        return df
