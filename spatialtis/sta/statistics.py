from collections import Counter
from typing import Optional, Sequence, Union

import numpy as np
import pandas as pd
from anndata import AnnData

from spatialtis.config import CONFIG
from spatialtis.utils import df2adata_uns, timer

Num = Union[int, float]


def type_counter(
    adata: AnnData,
    groupby: Union[Sequence, str, None] = None,
    type_key: Optional[str] = None,
    selected_types: Optional[Sequence] = None,
) -> pd.DataFrame:
    """(private) To count how many type of value

    Args:
        adata: anndata object to perform analysis
        groupby: how your experiments grouped, (Default: read from spatialtis.CONFIG.EXP_OBS)
        type_key: the key name of cell type in anndata.obs (Default: read from spatialtis.CONFIG.CELL_TYPE_KEY)
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

    if type_key is None:
        type_key = CONFIG.CELL_TYPE_KEY

    df = adata.obs[groupby + [type_key]]

    if selected_types is not None:
        df = df[df[type_key].isin(selected_types)]
    # the order of type will follow the order in df
    groups = df.groupby(groupby)
    types = pd.unique(df[type_key])
    matrix = list()
    mindex = list()
    for i, (n, g) in enumerate(groups):
        c = Counter(g[type_key])
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


@timer(prefix="Running cell components")
def cell_components(
    adata: AnnData,
    groupby: Union[Sequence, str, None] = None,
    type_key: Optional[str] = None,
    export: bool = True,
    export_key: Optional[str] = None,
    return_df: bool = False,
):
    """the proportion of different type of cells

    Args:
        adata: anndata object to perform analysis
        groupby: how your experiments grouped, (Default: read from spatialtis.CONFIG.EXP_OBS)
        type_key: the key name of cell type in anndata.obs (Default: read from spatialtis.CONFIG.CELL_TYPE_KEY)
        export: whether to export to anndata.uns field
        export_key: the key name that used to record the results in anndata.uns field
        return_df: whether to return an pandas.DataFrame

    Return:
        pandas.DataFrame

    """
    if export_key is None:
        export_key = CONFIG.cell_components_key
    else:
        CONFIG.cell_components_key = export_key

    counter = type_counter(adata, groupby, type_key)

    if export:
        df2adata_uns(counter, adata, export_key)

    if return_df:
        return counter


@timer(prefix="Running cell co-occurrence")
def cell_co_occurrence(
    adata: AnnData,
    groupby: Union[Sequence, str, None] = None,
    type_key: Optional[str] = None,
    export: bool = True,
    export_key: Optional[str] = None,
    threshold: int = 50,
    return_df: bool = False,
):
    """The probability of two type of cells occur simultaneously

        Args:
            adata: anndata object to perform analysis
            groupby: how your experiments grouped, (Default: read from spatialtis.CONFIG.EXP_OBS)
            type_key: the key name of cell type in anndata.obs (Default: read from spatialtis.CONFIG.CELL_TYPE_KEY)
            export: whether to export to anndata.uns field
            export_key: the key name that used to record the results in anndata.uns field
            threshold: this value determines the presence/absence of a cell type in ROI
            return_df: whether to return an pandas.DataFrame

        Return:
            pandas.DataFrame

    """
    if export_key is None:
        export_key = CONFIG.cell_co_occurrence_key
    else:
        CONFIG.cell_co_occurrence_key = export_key

    counter = type_counter(adata, groupby, type_key)

    occurrence = counter.gt(threshold)

    if export:
        df2adata_uns(occurrence, adata, export_key)

    if return_df:
        return occurrence


@timer(prefix="Running cell density")
def cell_density(
    adata: AnnData,
    size: Union[Sequence[Sequence[Num]], Sequence[Num]],
    ratio: float = 1.0,
    groupby: Union[Sequence, str, None] = None,
    type_key: Optional[str] = None,
    export: bool = True,
    export_key: Optional[str] = None,
    return_df: bool = False,
):
    """Calculating cell density in each ROI

    Args:
        adata: anndata object to perform analysis
        size: the size of each ROI, if some ROI is different, please specific each ROI size in a list, the order must
        follow the index order in anndata object.
        ratio: the ratio between pixel and real size, default is 1.0
        groupby: how your experiments grouped, (Default: read from spatialtis.CONFIG.EXP_OBS)
        type_key: the key name of cell type in anndata.obs (Default: read from spatialtis.CONFIG.CELL_TYPE_KEY)
        export: whether to export to anndata.uns field
        export_key: the key name that used to record the results in anndata.uns field (Default: "cell_density")
        return_df: whether to return an pandas.DataFrame

    Return:
            pandas.DataFrame

    """
    if export_key is None:
        export_key = CONFIG.cell_density_key
    else:
        CONFIG.cell_density_key = export_key

    counter = type_counter(adata, groupby, type_key)
    if isinstance(size[0], (int, float)):
        area = (size[0] * ratio) * (size[1] * ratio)
        results = counter.div(area)

    elif isinstance(size[0], Sequence):
        area = pd.Series({"area": [s[0] * s[1] * ratio * ratio for s in size]})
        results = counter.div(area["area"], axis=0)
    else:
        raise ValueError("Unrecognized size input")

    if export:
        df2adata_uns(results, adata, export_key)

    if return_df:
        return results


@timer(prefix="Running cell morphology")
def cell_morphology(
    adata: AnnData,
    metrics_key: str = "eccentricity",
    groupby: Union[Sequence, str, None] = None,
    type_key: Optional[str] = None,
    export: bool = True,
    export_key: Optional[str] = None,
    return_df: bool = False,
):
    """Cell morphology variation between different groups

        Args:
            adata: anndata object to perform analysis
            metrics_key: which key was used to measure cell morphology
            (Default: "eccentricity", Options: "eccentricity", "area")
            groupby: how your experiments grouped, (Default: read from spatialtis.CONFIG.EXP_OBS)
            type_key: the key name of cell type in anndata.obs (Default: read from spatialtis.CONFIG.CELL_TYPE_KEY)
            export: whether to export to anndata.uns field
            export_key: the key name that used to record the results in anndata.uns field (Default: "cell_morphology")
            return_df: whether to return an pandas.DataFrame

        Return:
                pandas.DataFrame

    """

    if groupby is None:
        groupby = CONFIG.EXP_OBS
    else:
        groupby = list(groupby)

    if type_key is None:
        type_key = CONFIG.CELL_TYPE_KEY

    if export_key is None:
        export_key = CONFIG.cell_morphology_key
    else:
        CONFIG.cell_morphology_key = export_key

    key = groupby + [type_key, metrics_key]
    df = adata.obs[key]

    df.index = pd.MultiIndex.from_frame(df[groupby + [type_key]])
    df = df.drop(groupby + [type_key], axis=1)

    if export:
        df2adata_uns(df, adata, export_key)

    if return_df:
        return df
