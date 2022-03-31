import warnings
from typing import Optional

import numpy as np
import pandas as pd
from anndata import AnnData
from spatialtis_core import comb_bootstrap

from spatialtis.abc import AnalysisBase
from spatialtis.typing import Array
from spatialtis.utils import doc
from spatialtis.utils.io import read_neighbors


@doc
def spatial_enrichment(data: AnnData,
                       threshold: Optional[float] = None,
                       layer_key: Optional[str] = None,
                       selected_markers: Optional[Array] = None,
                       resample: int = 500,
                       pval: float = 0.01,
                        export_key: str = "spatial_enrichment",
                       **kwargs, ):
    """`Profiling markers spatial enrichment <about/implementation.html#profiling-of-markers-co-expression>`_
    using permutation test

    Similar to neighborhood analysis which tells you the relationship between different type of cells.
    This analysis tells you the spatial relationship between markers.

    Args:
        data: {adata}
        threshold: The expression level to determine whether a marker is positive
        layer_key: {layer_key}
        selected_markers: {selected_markers}
        resample: Number of times to perform resample
        pval: {pval}
        export_key: {export_key}
        **kwargs: {analysis_kwargs}

    .. seealso:: :class:`spatialtis.cell_interaction`

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
    for roi_name, roi_data, mks, exp in ab.roi_exp_iter(
            selected_markers=markers,
            layer_key=layer_key,
            dtype=np.bool,
            desc="Spatial enrichment",
    ):
        neighbors = read_neighbors(roi_data, ab.neighbors_key)
        labels = roi_data[ab.cell_id_key]
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
