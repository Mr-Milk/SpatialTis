from collections import Counter

import numpy as np
import pandas as pd
from anndata import AnnData
from spatialtis_core import somde as smode_sv

from spatialtis.abc import AnalysisBase
from spatialtis.utils import doc


@doc
def somde(data: AnnData,
          k: int = 20,
          alpha: float = 0.5,
          epoch: int = 100,
          pval: float = 0.05,
          qval: float = 0.05,
          **kwargs, ):
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
    ab = AnalysisBase(data, display_name="SOMDE", **kwargs)
    track_ix = []
    results_data = []
    for roi_name, roi_data, markers, exp, points in ab.roi_exp_iter_with_points(
            desc="Spatial variable genes: SOMDE"
    ):
        sv_genes = smode_sv(
            pd.DataFrame(exp, index=markers, dtype=np.float32).fillna(0.0),
            np.array(points, dtype=np.float32),
            k=k, alpha=alpha, epoch=epoch, pval=pval, qval=qval
        )
        results_data.append(sv_genes)
        track_ix.append(roi_name)

    # a dict store all the markers
    markers_dict = {k: 0 for k in ab.markers}
    # unpack and merge to ensure every counter has the same markers
    ab.result = pd.DataFrame(
        data=[{**markers_dict, **Counter(i)} for i in results_data],
        index=pd.MultiIndex.from_tuples(track_ix, names=ab.exp_obs),
    )[ab.markers]
