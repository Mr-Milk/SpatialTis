from __future__ import annotations

import numpy as np
import pandas as pd
import scipy
from anndata import AnnData
from ast import literal_eval
from scipy.stats import mannwhitneyu
from spatialtis_core import neighbor_components
from typing import Dict, List

from spatialtis.abc import AnalysisBase
from spatialtis.utils import doc, read_exp


@doc
def NCD_marker(data: AnnData,
               selected_markers: List[str] | np.ndarray = None,
               importance_cutoff: float = 0.5,
               layer_key: str = None,
               tree_kwargs: Dict = None,
               test_method: str = "mannwhitneyu",
               pval: float = 0.01,
               export_key: str = "ncd_marker",
               **kwargs, ):
    """Identify neighbor cells dependent marker

    This method tells you the dependency between markers and its neighbor cell type.
    The dependency is calculated by building a gradiant boosting tree (in here lightgbm) to determine
    the feature importance. A statistic test and fold change will be calculated for importance markers and its
    neighbor cells, the fold change is between marker with cell type at / not at the neighborhood.

    Parameters
    ----------
    data : {adata}
    importance_cutoff : float, default: 0.5
        Threshold to determine the feature markers.
    selected_markers : {selected_markers}
    layer_key : {layer_key}
    tree_kwargs : dict
        The keyword arguments that pass to the boosting tree class, (Default: n_jobs=-1, random_state=0).
    test_method : str, default: 'mannwhitneyu'
        which test method to use, anything from
        `scipy.stats <https://docs.scipy.org/doc/scipy/reference/stats.html>`_.
    pval : {pval}
    export_key : {export_key}
    **kwargs : {analysis_kwargs}

    """

    try:
        from lightgbm import LGBMRegressor
    except ImportError:
        raise ImportError("lightgbm is not installed, please try `pip install lightgbm`.")
    ab = AnalysisBase(data,
                      display_name="NCD Markers",
                      export_key=export_key,
                      **kwargs)
    ab.check_neighbors()
    ab.check_cell_type()

    tree_kwargs_ = {"n_jobs": -1, "random_state": 0, "importance_type": "gain"}
    if tree_kwargs is not None:
        for k, v in tree_kwargs.items():
            tree_kwargs_[k] = v

    markers = ab.selected_markers(selected_markers)
    markers_mask = ab.markers_col.isin(markers)

    neighbors = [literal_eval(n) for n in data.obsm[ab.neighbors_key]]
    labels = data.obs[ab.cell_id_key]
    cell_types = data.obs[ab.cell_type_key]
    col, comps = neighbor_components(
        neighbors, labels.tolist(), cell_types.tolist()
    )
    neigh_comp = pd.DataFrame(
        data=comps,
        columns=col,
        index=pd.MultiIndex.from_frame(
            data.obs[[ab.cell_type_key, ab.cell_id_key]],
            names=["type", "id"],
        ),
    )
    results_data = []
    # For markers in different cell types
    with np.errstate(divide="ignore"):
        for t, x in neigh_comp.groupby(level=["type"]):
            exp_ix = x.index.to_frame()["id"]
            exp = read_exp(data[exp_ix, markers_mask], layer_key)
            for i, y in enumerate(exp):
                # copy it to prevent memory peak according to lightgbm
                reg = LGBMRegressor(**tree_kwargs_).fit(x, y.copy())
                weights = np.asarray(reg.feature_importances_)
                weights = weights / weights.sum()
                max_ix = np.argmax(weights)
                max_weight = weights[max_ix]
                max_type = col[max_ix]
                if max_weight > importance_cutoff:
                    nx = x.copy()
                    # add expression data to dataframe to allow cutting afterwards
                    nx["exp"] = y
                    # cells with max_type at neighbors
                    at_neighbor = (nx.iloc[:, max_ix] != 0)
                    at_neighbor_exp = nx[at_neighbor]["exp"].to_numpy()
                    non_at_neighbor_exp = nx[~at_neighbor]["exp"].to_numpy()
                    at_sum = at_neighbor_exp.sum()
                    non_at_sum = non_at_neighbor_exp.sum()
                    if (at_sum > 0) & (non_at_sum > 0):
                        test_result = getattr(scipy.stats, test_method).__call__(
                            at_neighbor_exp, non_at_neighbor_exp
                        )
                        pvalue = test_result.pvalue
                        if pvalue < pval:
                            at_mean = at_neighbor_exp.mean()
                            non_at_mean = non_at_neighbor_exp.mean()
                            log2_fc = np.log2(at_mean / non_at_mean)
                            results_data.append([t, markers[i], max_type,
                                                 max_weight, log2_fc, pvalue, ])
    ab.result = pd.DataFrame(
        data=results_data,
        columns=[
            "cell_type",
            "marker",
            "neighbor_type",
            "dependency",
            "log2_FC",
            "pval",
        ],
    )
