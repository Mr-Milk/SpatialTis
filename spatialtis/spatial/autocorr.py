from typing import Optional

import numpy as np
import pandas as pd
from anndata import AnnData
from spatialtis_core import spatial_autocorr as autocorr

from spatialtis.abc import AnalysisBase
from spatialtis.utils import read_neighbors, doc, options_guard


@doc
def spatial_autocorr(
        data: AnnData,
        method: str = "moran_i",
        pval: float = 0.05,
        two_tailed: bool = True,
        layer_key: Optional[str] = None,
        export_key: str = "spatial_autocorr",
        **kwargs,
):
    """Spatial auto-correlation for every markers

    This is used measure the correlation of marker expression with spatial locations.

    Moran's I is more for global spatial autocorrelation,
    Geary's C is more for local spatial autocorrelation

    Args:
        data: {data}
        method: "moran_i" or "geary_c" (Default: "moran_i")
        pval: {pval}
        two_tailed: Whether to use two tailed for p-value
        layer_key: {layer_key}
        export_key: {export_key}
        **kwargs: {analysis_kwargs}

    .. seealso:: :class:`spatialtis.somde`

    """
    method = options_guard(method, ['moran_i', 'geary_c'])
    ab = AnalysisBase(data,
                      method=method,
                      display_name="Spatial auto-correlation",
                      export_key=export_key,
                      **kwargs)
    track_ix = []
    results_data = []
    for roi_name, roi_data, markers, exp in ab.roi_exp_iter(
            layer_key=layer_key, desc=ab.display_name
    ):
        neighbors = read_neighbors(roi_data, ab.neighbors_key)
        labels = roi_data[ab.cell_id_key]
        results = autocorr(
            exp.astype(np.float64),
            neighbors,
            labels=labels,
            two_tailed=two_tailed,
            pval=pval,
            method=method,
        )
        markers = markers.to_numpy()
        results = np.hstack([markers.reshape(-1, 1), results])
        track_ix += [roi_name for _ in range(len(markers))]
        results_data.append(results)

    ab.result = pd.concat(
        [
            pd.DataFrame(data=track_ix, columns=ab.exp_obs),
            pd.DataFrame(
                data=np.concatenate(results_data),
                columns=["marker", "pattern", "index_value", "pval"],
            ),
        ],
        axis=1,
    )
