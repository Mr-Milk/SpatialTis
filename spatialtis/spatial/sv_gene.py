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
          export_key: str = "sv_gene",
          **kwargs, ):
    """This is a wrapper around somde

    Parameters
    ----------
    data : {adata}
    k : int, default: 20
        Number of SOM nodes.
    alpha : float, default: 0.5
        Parameters for generate pseudo gene expression.
    epoch : int, default: 100
        Number of epoch.
    qval : float, default: 0.05
        Threshold for qval.
    pval : float, default: 0.05
        Threshold for pval.
    export_key : {export_key}
    **kwargs : {analysis_kwargs}

    """
    ab = AnalysisBase(data, display_name="SOMDE", export_key=export_key, **kwargs)
    track_ix = []
    results_data = []
    for roi_name, markers, exp, points in ab.iter_roi(
            fields=['exp', 'centroid'], desc="Spatial variable genes: SOMDE"
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
