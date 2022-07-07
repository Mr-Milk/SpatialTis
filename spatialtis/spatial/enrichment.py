from __future__ import annotations

import numpy as np
import pandas as pd
import warnings
from anndata import AnnData
from spatialtis_core import comb_bootstrap
from typing import List

from spatialtis.abc import AnalysisBase
from spatialtis.utils import doc


@doc
def spatial_enrichment(data: AnnData,
                       threshold: float = None,
                       layer_key: str = None,
                       selected_markers: List[str] | np.ndarray = None,
                       resample: int = 500,
                       pval: float = 0.01,
                       export_key: str = "spatial_enrichment",
                       **kwargs, ):
    """`Profiling markers spatial enrichment <about/implementation.html#profiling-of-markers-co-expression>`_
    using permutation test

    Similar to neighborhood analysis which tells you the relationship between different type of cells.
    This analysis tells you the spatial relationship between markers.

    Parameters
    ----------
    data : {adata}
    threshold : float
        The expression level to determine whether a marker is positive
    layer_key : {layer_key}
    selected_markers : {selected_markers}
    resample : int, default: 500
        Number of times to perform resample
    pval : {pval}
    export_key : {export_key}
    **kwargs : {analysis_kwargs}

    See Also
    --------
    :class:`spatialtis.cell_interaction`

    """

    ab = AnalysisBase(data,
                      display_name="Spatial enrichment",
                      export_key=export_key,
                      **kwargs)
    ab.check_neighbors()

    if (threshold is not None) & (layer_key is None):
        layer_key = f"gt_{threshold}"
        data.layers[layer_key] = (data.X.copy() >= threshold).astype(bool)
    elif (threshold is not None) & (layer_key is not None):
        warnings.warn(
            "You specific both threshold and layers_key, "
            "using user defined layers_key"
        )
    else:
        layer_key = f"mean_cut"
        data.layers[layer_key] = (data.X.copy() >= data.X.mean(axis=0)).astype(
            bool
        )

    markers = ab.selected_markers(selected_markers)

    results_data = []
    for roi_name, mks, exp, labels, neighbors in ab.iter_roi(
            fields=["exp", "neighbors"],
            selected_markers=markers,
            layer_key=layer_key,
            dtype=np.bool,
    ):
        result = comb_bootstrap(
            exp,
            mks,
            neighbors,
            labels,
            pval=pval,
            times=resample,
        )
        for pairs in result:
            results_data.append([*roi_name, *pairs])

    df = pd.DataFrame(
        data=results_data, columns=ab.exp_obs + ["marker1", "marker2", "value"]
    )
    df = df.pivot_table(
        values="value", index=ab.exp_obs, columns=["marker1", "marker2"]
    )
    ab.result = df
