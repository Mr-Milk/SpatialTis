from collections import Counter
from typing import Optional

import pandas as pd
from anndata import AnnData
from spatialtis_core import somde as smode_sv

from spatialtis.abc import AnalysisBase
from spatialtis.utils import read_points


class somde(AnalysisBase):
    def __init__(self,
                 data: AnnData,
                 marker_key: Optional[str] = None,
                 k: int = 20,
                 alpha: float = 0.5,
                 epoch: int = 100,
                 pval: float = 0.05,
                 qval: float = 0.05,
                 ):
        super().__init__(self, data)
        track_ix = []
        results_data = []
        for roi_name, roi_data, markers, exp in self.roi_exp_iter(desc="Spatial variable genes: SOMDE"):
            points = read_points(roi_data, self.centroid_key)
            sv_genes = smode_sv(exp, points, k=k, alpha=alpha, epoch=epoch, pval=pval, qval=qval)
            results_data.append(sv_genes)
            track_ix.append(roi_name)

        # a dict store all the markers
        markers_dict = {k:0 for k in self.markers}
        # unpack and merge to ensure every counter has the same markers
        self.results = pd.DataFrame(data=[{**markers_dict, **Counter(i)} for i in results_data],
                                    index=pd.MultiIndex.from_tuples(track_ix))[self.markers]
