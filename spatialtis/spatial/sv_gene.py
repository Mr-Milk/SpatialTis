from collections import Counter
from typing import Optional

import numpy as np
import pandas as pd
from anndata import AnnData
from spatialtis_core import somde as smode_sv

from spatialtis.abc import AnalysisBase
from spatialtis.utils import read_points, doc


@doc
class somde(AnalysisBase):
    """This is a wrapper around somde

    Args:
        data: {adata}
        k: Number of SOM nodes
        alpha: Parameters for generate pseudo gene expression
        epoch: Number of epoch
        qval: Threshold for qval
        pval: Threshold for pval
        **kwargs: {analysis_kwargs}

    """
    def __init__(
        self,
        data: AnnData,
        k: int = 20,
        alpha: float = 0.5,
        epoch: int = 100,
        pval: float = 0.05,
        qval: float = 0.05,
        **kwargs,
    ):
        self.display_name = "SOMDE"
        super().__init__(data, **kwargs)
        track_ix = []
        results_data = []
        for roi_name, roi_data, markers, exp in self.roi_exp_iter(
            desc="Spatial variable genes: SOMDE"
        ):
            points = read_points(roi_data, self.centroid_key)
            sv_genes = smode_sv(
                pd.DataFrame(exp, index=markers, dtype=np.float32).fillna(0.0),
                np.array(points, dtype=np.float32),
                k=k, alpha=alpha, epoch=epoch, pval=pval, qval=qval
            )
            results_data.append(sv_genes)
            track_ix.append(roi_name)

        # a dict store all the markers
        markers_dict = {k: 0 for k in self.markers}
        # unpack and merge to ensure every counter has the same markers
        self.result = pd.DataFrame(
            data=[{**markers_dict, **Counter(i)} for i in results_data],
            index=pd.MultiIndex.from_tuples(track_ix, names=self.exp_obs),
        )[self.markers]
