import warnings
from typing import Optional, Sequence

import numpy as np
import pandas as pd
from tqdm import tqdm

from spatialtis.config import CONFIG
from spatialtis.spatial.neighbors import Neighbors
from spatialtis.spatial.utils import check_neighbors
from spatialtis.utils import df2adata_uns, get_default_params, reuse_docstring, timer


@timer(prefix="Running spatial enrichment analysis")
@get_default_params
@reuse_docstring()
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
    """`Profiling markers co-expression <about/implementation.html#profiling-of-markers-co-expression>`_ using permutation test

    Similar to neighborhood analysis which tells you the relationship between different type of cells.
    This analysis tells you the spatial relationship between markers.

    This method is implemented in Rust, it executes in parallel automatically.

    Args:
        n: {n}
        threshold: The expression level to determine whether a marker is positive
        layers_key: {layers_key}
        selected_markers: {selected_markers}
        marker_key: {marker_key}
        resample: Number of times to perform resample
        pval: {pval}
        export: {export}
        export_key: {export_key}
        return_df: {return_df}

    .. seealso:: `neighborhood_analysis <#spatialtis.plotting.neighborhood_analysis>`_

    """

    try:
        import neighborhood_analysis as na
    except ImportError:
        raise ImportError("Package not found, try `pip install neighborhood_analysis`.")

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
            warnings.warn(
                "You specific both threshold and layers_key, using user defined layers_key"
            )
        CONFIG.spatial_enrichment_analysis_layers_key = layers_key
    elif threshold is not None:
        layers_key = CONFIG.spatial_enrichment_analysis_layers_key
        data.layers[layers_key] = data.X >= threshold

    if selected_markers is not None:
        if len(selected_markers) > 1:
            data_t = data.T
            data = data_t[data_t.obs[marker_key].isin(selected_markers)].copy()
            data = data.T
        else:
            raise ValueError("You need at least two markers for `selected_markers`.")

    markers = data.var[marker_key]
    results = {}

    for name, roi in tqdm(
        data.obs.groupby(n.expobs), **CONFIG.tqdm(desc="Spatial enrichment analysis"),
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
