import warnings
from typing import Optional, Sequence

import numpy as np
import pandas as pd
from tqdm import tqdm

from spatialtis.config import CONFIG
from spatialtis.utils import df2adata_uns, timer

from ._neighbors import Neighbors
from ._util import check_neighbors


@timer(prefix="Running neighborhood analysis")
def neighborhood_analysis(
    n: Neighbors,
    method: str = "pval",
    resample: int = 500,
    pval: float = 0.01,
    order: bool = True,
    export: bool = True,
    export_key: Optional[str] = None,
    return_df: bool = False,
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

    .. seealso:: `spatial_enrichment_analysis <#spatialtis.plotting.spatial_enrichment_analysis>`_


    """

    if export_key is None:
        export_key = CONFIG.neighborhood_analysis_key
    else:
        CONFIG.neighborhood_analysis_key = export_key

    try:
        import neighborhood_analysis as na
    except ImportError:
        raise ImportError("Package not found, try `pip install neighborhood_analysis`.")

    check_neighbors(n)
    types = n.unitypes
    cc = na.CellCombs(types, order)

    results = {}
    for name, value in tqdm(
        n.neighbors.items(), **CONFIG.tqdm(desc="neighborhood analysis"),
    ):
        result = cc.bootstrap(
            n.types[name], value, resample, pval, method, ignore_self=True
        )
        result = {tuple(k): v for (k, v) in result}
        results[name] = result

    df = pd.DataFrame(results)
    df.index = pd.MultiIndex.from_tuples(df.index, names=("Cell type1", "Cell type2"))
    df.rename_axis(columns=n.expobs, inplace=True)

    if method == "pval":
        df = df.T.astype(int)
    else:
        df.replace(np.inf, 10.0, inplace=True)
        df.replace(-np.inf, -10.0, inplace=True)
        df = df.T

    if export:
        df2adata_uns(
            df,
            n.adata,
            export_key,
            params={"order": order, "method": method, "pval": pval},
        )

    if return_df:
        return df


# @timer(prefix="Running spatial enrichment analysis")
def spatial_enrichment_analysis(
    n: Neighbors,
    threshold: Optional[float] = None,
    layers_key: Optional[str] = None,
    selected_markers: Optional[Sequence] = None,
    marker_key: Optional[str] = None,
    resample: int = 500,
    pval: float = 0.01,
    export: bool = True,
    export_key: Optional[str] = None,
    return_df: bool = False,
):
    """Profiling Markers co-expression

    Similar to neighborhood analysis which tells you the relationship between different type of cells.
    This method tells you the spatial relationship between markers, purposed in MIBI's paper.

    Args:
        n: A spatialtis.Neighbors object, neighbors are already computed
        threshold: the number to determine whether a marker is positive
        layers_key: the key to anndata.layers
        selected_markers: the markers to perform analysis on
        marker_key: the key of markers in anndata.var (Default: spatialtis.CONFIG.MARKER_KEY)
        resample: perform resample for how many times
        pval: threshold for p-value
        export: whether to export the result to anndata.uns
        export_key: the key used to export
        return_df: whether to return the result

    .. seealso:: `neighborhood_analysis <#spatialtis.plotting.neighborhood_analysis>`_

    """

    try:
        import neighborhood_analysis as na
    except ImportError:
        raise ImportError("Package not found, try `pip install neighborhood_analysis`.")

    check_neighbors(n)
    data = n.adata

    if marker_key is None:
        marker_key = CONFIG.MARKER_KEY

    if export_key is None:
        export_key = CONFIG.spatial_enrichment_analysis_key
    else:
        CONFIG.spatial_enrichment_analysis_key = export_key

    if (threshold is None) & (layers_key is None):
        raise ValueError("Either specific a threshold or a layers key.")
    elif layers_key is not None:
        if threshold is not None:
            warnings.warn(
                "You specific both threshold and layers_key, using user defined layers_key"
            )
        CONFIG.spatial_enrichment_analysis_layers_key = layers_key
    elif threshold is not None:
        layers_key = CONFIG.spatial_enrichment_analysis_layers_key
        data.layers[layers_key] = data.X >= threshold

    if selected_markers is not None:
        data = data[data.var[marker_key].isin(selected_markers)]

    markers = data.var[marker_key]
    results = {}

    for name, roi in tqdm(
        data.obs.groupby(n.expobs), **CONFIG.tqdm(desc="spatial enrichment analysis"),
    ):
        neighbors = n.neighbors[name]
        matrix = data[roi.index].layers[layers_key]
        result = {}

        for ix, x in enumerate(markers):
            x_status = [bool(i) for i in matrix[:, ix]]
            for iy, y in enumerate(markers):
                y_status = [bool(i) for i in matrix[:, iy]]
                z = na.comb_bootstrap(
                    x_status, y_status, neighbors, times=resample, ignore_self=False
                )
                result[(x, y)] = z
        results[name] = result

    df = pd.DataFrame(results)
    df.index = pd.MultiIndex.from_tuples(df.index, names=("Marker1", "Marker2"))
    df.rename_axis(columns=n.expobs, inplace=True)
    df.replace(np.inf, 10.0, inplace=True)
    df.replace(-np.inf, -10.0, inplace=True)
    df.fillna(0, inplace=True)
    df = df.T

    if export:
        df2adata_uns(df, n.adata, export_key, params={"pval": pval})

    if return_df:
        return df
