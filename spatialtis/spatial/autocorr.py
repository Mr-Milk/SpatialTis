from typing import Optional

import numpy as np
import pandas as pd
from anndata import AnnData
from spatialtis_core import spatial_autocorr as autocorr

from spatialtis.abc import AnalysisBase
from spatialtis.utils import read_neighbors, doc


@doc
class spatial_autocorr(AnalysisBase):
    """Spatial auto-correlation for every markers

    Args:
        data: {data}
        marker_key: {marker_key}
        layer_key: {layer_key}
        two_tailed: If using
        method:
    """
    def __init__(
        self,
        data: AnnData,
        marker_key: Optional[str] = None,
        layer_key: Optional[str] = None,
        two_tailed: bool = True,
        method: str = "moran_i",
    ):
        super().__init__(data, marker_key=marker_key)
        track_ix = []
        results_data = []
        for roi_name, roi_data, markers, exp in self.roi_exp_iter(
            layer_key=layer_key, desc="Spatial autocorrelation"
        ):
            neighbors = read_neighbors(roi_data, self.neighbors_key)
            labels = roi_data[self.cell_id_key]
            results = autocorr(
                exp.astype(np.float64),
                neighbors,
                labels=labels,
                two_tailed=two_tailed,
                method=method,
            )
            markers = markers.to_numpy()
            results = np.hstack([markers.reshape(-1, 1), results])
            track_ix += [roi_name for _ in range(len(markers))]
            results_data.append(results)

        self.result = pd.concat(
            [
                pd.DataFrame(data=track_ix, columns=self.exp_obs),
                pd.DataFrame(
                    data=np.concatenate(results_data),
                    columns=["marker", "pattern", "index_value", "pval"],
                ),
            ],
            axis=1,
        )
