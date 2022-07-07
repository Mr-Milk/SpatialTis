from __future__ import annotations

import numpy as np
import pandas as pd
from anndata import AnnData
from ast import literal_eval
from scipy.stats import spearmanr
from typing import Dict, List

from spatialtis.abc import AnalysisBase
from spatialtis.utils import doc, pbar_iter, read_exp


@doc
def NMD_marker(data: AnnData,
               pval: float = 0.01,
               selected_markers: List[str] | np.ndarray = None,
               importance_cutoff: float = 0.5,
               layer_key: str = None,
               tree_kwargs: Dict = None,
               export_key: str = "nmd_marker",
               **kwargs, ):
    """Identify neighbor markers dependent marker

    The neighborhood is treated as a single cell.

    Parameters
    ----------
    data : {adata}
    importance_cutoff : float
        Standard deviation, threshold to filter out markers that are not variant enough.
    pval : {pval}
    selected_markers : {selected_markers}
    layer_key : {layers_key}
    tree_kwargs : dict
        The keyword arguments that pass to the boosting tree class, (Default: n_jobs=-1, random_state=0).
    export_key : {export_key}
    **kwargs : {analysis_kwargs}

    """
    try:
        from lightgbm import LGBMRegressor
    except ImportError:
        raise ImportError("lightgbm is not installed, please try `pip install lightgbm`.")
    ab = AnalysisBase(data, display_name="NMD marker", export_key=export_key, **kwargs)
    ab.check_neighbors()

    tree_kwargs_ = {"n_jobs": -1, "random_state": 0, "importance_type": "gain"}
    if tree_kwargs is not None:
        for k, v in tree_kwargs.items():
            tree_kwargs_[k] = v

    markers = ab.selected_markers(selected_markers)
    markers_mask = ab.markers_col.isin(markers)
    neighbors = [literal_eval(n) for n in data.obsm[ab.neighbors_key]]
    cent_exp = read_exp(data[:, markers_mask], layer_key)
    # treat the neighbors as single cell
    # sum the expression
    neigh_exp = np.asarray([read_exp(data[n, markers_mask], layer_key).sum(1) for n in neighbors])
    results_data = []
    for i, y in enumerate(pbar_iter(cent_exp, desc="Neighbor-dependent markers", )):
        reg = LGBMRegressor(**tree_kwargs_).fit(neigh_exp, y)
        weights = np.asarray(reg.feature_importances_)
        ws = weights.sum()
        if ws != 0:
            weights = weights / weights.sum()
            max_ix = np.argmax(weights)
            max_weight = weights[max_ix]
            if max_weight > importance_cutoff:
                r, pvalue = spearmanr(y, neigh_exp[:, max_ix])
                if pvalue < pval:
                    results_data.append(
                        [markers[i], markers[max_ix], max_weight, r, pvalue]
                    )

    ab.result = pd.DataFrame(
        data=results_data,
        columns=["marker", "neighbor_marker", "dependency", "corr", "pval"],
    )
