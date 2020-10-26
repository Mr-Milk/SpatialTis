from ast import literal_eval
from collections import Counter
from itertools import combinations_with_replacement
from typing import Optional, Sequence, Union

import numpy as np
import pandas as pd
from anndata import AnnData
from scipy.stats import chisquare
from shapely.geometry import MultiPoint

from spatialtis.config import CONFIG
from spatialtis.utils import df2adata_uns, get_default_params, reuse_docstring, timer

Num = Union[int, float]


@reuse_docstring()
def type_counter(
    adata: AnnData,
    groupby: Union[Sequence, str, None] = None,
    type_key: Optional[str] = None,
) -> pd.DataFrame:
    """(private) To count how many types of value

    Args:
        adata: {adata}
        groupby: {groupby}
        type_key: {type_key}

    """

    df = adata.obs[groupby + [type_key]]

    # the order of types will follow the order in df
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
@get_default_params
@reuse_docstring()
def cell_components(
    adata: AnnData,
    groupby: Union[Sequence, str, None] = None,
    type_key: Optional[str] = None,
    export: bool = True,
    export_key: Optional[str] = None,
    return_df: bool = False,
):
    """Count the proportion of each types of cells in each group

    Args:
        adata: {adata}
        groupby: {groupby}
        type_key: {type_key}
        export: {export}
        export_key: {export_key}
        return_df: {return_df}

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
@get_default_params
@reuse_docstring()
def cell_co_occurrence(
    adata: AnnData,
    groupby: Union[Sequence, str, None] = None,
    threshold: int = 50,
    pval: float = 0.01,
    type_key: Optional[str] = None,
    export: bool = True,
    export_key: Optional[str] = None,
    return_df: bool = False,
):
    """The likelihood of two type of cells occur simultaneously in a ROI

    Using chi-square test to determine the significance of co-occurrence.

    Args:
        adata: {adata}
        groupby: {groupby}
        threshold: The threshold value (number of cells) determines the presence/absence of a cell type in ROI
        pval: {pval}
        type_key: {type_key}
        export: {export}
        export_key: {export_key}
        return_df: {return_df}

    """
    if export_key is None:
        export_key = CONFIG.cell_co_occurrence_key
    else:
        CONFIG.cell_co_occurrence_key = export_key

    counter = type_counter(adata, groupby, type_key)

    df = counter.gt(threshold)

    groups = df.astype(int).groupby(level=groupby)
    X = []
    for n, g in groups:
        c = (g.sum() / len(g)).fillna(0)
        data = []
        for comb in combinations_with_replacement(c, 2):
            if (comb[0] == 0) | (comb[1] == 0):
                p = 0
            else:
                p = chisquare(comb).pvalue
            data.append(p)
        X.append(data)
    pdf = (pd.DataFrame(X) > (1 - pval)).astype(int)
    pdf.index = pd.MultiIndex.from_frame(
        df.index.to_frame(index=False)[groupby].drop_duplicates()
    )
    pdf.columns = pd.MultiIndex.from_arrays(
        np.asarray([i for i in combinations_with_replacement(df.columns, 2)]).T,
        names=["Cell type1", "Cell type2"],
    )

    if export:
        df2adata_uns(pdf, adata, export_key)

    if return_df:
        return pdf


@timer(prefix="Running cell density")
@get_default_params
@reuse_docstring()
def cell_density(
    adata: AnnData,
    groupby: Union[Sequence, str, None] = None,
    size: Union[Sequence[Sequence[Num]], Sequence[Num], None] = None,
    ratio: float = 1.0,
    type_key: Optional[str] = None,
    centroid_key: Optional[str] = None,
    export: bool = True,
    export_key: Optional[str] = None,
    return_df: bool = False,
):
    """Calculating cell density in each ROI

    Args:
        adata: {adata}
        groupby: {groupby}
        size: The size of each ROI.
              If some ROI is different, please specific each ROI size in a list,
              the order must follow the index order in AnnData object;
              If None, it will automatically computed the area of every ROI;
        ratio: The ratio between 1 pixel and real size, default is 1.0;
              For example, if your density unit is n cells/mm^2, your resolution is 1μm,
              then you should set the ratio as 0.001, 1 pixels represent 0.001mm length.
        type_key: {type_key}
        centroid_key: {centroid_key}
        export: {export}
        export_key: {export_key}
        return_df: {return_df}

    """
    if export_key is None:
        export_key = CONFIG.cell_density_key
    else:
        CONFIG.cell_density_key = export_key

    counter = type_counter(adata, groupby, type_key)

    if size is None:
        groups = adata.obs[groupby + [centroid_key]].groupby(groupby)
        area = [
            MultiPoint(
                [literal_eval(c) for c in g[centroid_key].tolist()]
            ).convex_hull.area
            for _, g in groups
        ]
        area = np.asarray(area) * ratio * ratio
        results = counter.div(area, axis=0)
    elif isinstance(size[0], (int, float)):
        area = (size[0] * ratio) * (size[1] * ratio)
        results = counter.div(area)

    elif isinstance(size[0], Sequence):
        area = pd.Series({"area": [s[0] * s[1] * ratio * ratio for s in size]})
        results = counter.div(area["area"], axis=0)
    else:
        raise TypeError("Type of parameter `size` incorrect.")

    if export:
        df2adata_uns(results, adata, export_key)

    if return_df:
        return results


@timer(prefix="Running cell morphology")
@get_default_params
@reuse_docstring()
def cell_morphology(
    adata: AnnData,
    groupby: Union[Sequence, str, None] = None,
    metric_key: Optional[str] = None,
    type_key: Optional[str] = None,
    export: bool = True,
    export_key: Optional[str] = None,
    return_df: bool = False,
):
    """Cell morphology variation between different groups

    Args:
        adata: {adata}
        groupby: {groupby}
        metric_key: {metric_key}
        type_key: {type_key}
        export: {export}
        export_key: {export_key}
        return_df: {return_df}

    """
    if export_key is None:
        export_key = CONFIG.cell_morphology_key
    else:
        CONFIG.cell_morphology_key = export_key

    key = groupby + [type_key, metric_key]
    df = adata.obs.loc[:, key]
    df = df.reset_index()
    df = df.rename(columns={type_key: "type", metric_key: "value", "index": "id"})
    df = df.set_index(groupby + ["type", "id"], drop=True)

    if export:
        df2adata_uns(df, adata, export_key)

    if return_df:
        return df
