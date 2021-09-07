from typing import Optional

import numpy as np
import pandas as pd
from anndata import AnnData
from spatialtis_core import spatial_autocorr as autocorr

from spatialtis.abc import AnalysisBase
from spatialtis.utils import read_neighbors


class spatial_autocorr(AnalysisBase):
    def __init__(self,
                 data: AnnData,
                 marker_key: Optional[str] = None,
                 layer_key: Optional[str] = None,
                 two_tailed=True,
                 method="moran_i"):
        super().__init__(data, marker_key=marker_key)
        track_ix = []
        markers_list = []
        results_data = []
        for roi_name, roi_data, markers, exp in self.roi_exp_iter(layer_key=layer_key,
                                                                  desc="Spatial autocorrelation"):
            neighbors = read_neighbors(roi_data, self.neighbors_key)
            labels = roi_data[self.cell_id_key]
            results = autocorr(exp, neighbors, labels=labels, two_tailed=two_tailed, method=method)
            track_ix.append(roi_name)
            markers_list.append(markers)
            results_data.append(results)

        self.result = pd.DataFrame(data=np.hstack([np.array(markers_list), np.array(results_data)]),
                                   columns=['marker', 'pattern', 'index_value', 'pval'],
                                   index=pd.MultiIndex.from_tuples(track_ix))