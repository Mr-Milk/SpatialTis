import numpy as np
import pandas as pd
from anndata import AnnData
from spatialtis_core import spatial_autocorr as autocorr
from typing import Literal

from spatialtis.abc import AnalysisBase
from spatialtis.utils import doc, options_guard


@doc
def spatial_autocorr(
        data: AnnData,
        method: Literal["moran_i", "geary_c"] = "moran_i",
        pval: float = 0.05,
        two_tailed: bool = True,
        layer_key: str = None,
        export_key: str = "spatial_autocorr",
        **kwargs,
):
    """Spatial auto-correlation for every marker.

    Measure the correlation of marker expression with spatial locations.

    - **Moran's I** is more for global spatial autocorrelation.
    - **Geary's C** is more for local spatial autocorrelation.

    Parameters
    ----------
    data : {adata}
    method : {'moran_i', 'geary_c'}, default: 'moran_i'
    pval : {pval}
    two_tailed : bool, default: True
        If True, p-value is two-tailed.
    layer_key : {layer_key}
    export_key : {export_key}
    **kwargs : {analysis_kwargs}

    See Also
    --------
    :class:`spatialtis.somde`

    """
    method = options_guard(method, ['moran_i', 'geary_c'])
    ab = AnalysisBase(data,
                      method=method,
                      display_name="Spatial auto-correlation",
                      export_key=export_key,
                      **kwargs)
    track_ix = []
    results_data = []
    for roi_name, labels, neighbors, markers, exp in ab.iter_roi(
            fields=['neighbors', 'exp'], layer_key=layer_key
    ):
        results = autocorr(
            exp.astype(np.float64),
            neighbors,
            labels=labels,
            two_tailed=two_tailed,
            pval=pval,
            method=method,
        )
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
