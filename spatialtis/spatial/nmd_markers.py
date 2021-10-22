from typing import Dict, Optional

import numpy as np
import pandas as pd
from anndata import AnnData
from scipy.stats import spearmanr

from spatialtis.abc import AnalysisBase
from spatialtis.utils import NeighborsNotFoundError
from spatialtis.typing import Number
from spatialtis.utils import doc, pbar_iter, read_exp, read_neighbors


@doc
class NMD_marker(AnalysisBase):
    """Identify neighbor markers dependent marker

    The neighborhood is treated as a single cell.

    Args:
        data: {adata}
        exp_std_cutoff: Standard deviation, threshold to filter out markers that are not variant enough
        pval: {pval}
        selected_markers: {selected_markers}
        layers_key: {layers_key}
        tree_kwargs: {tree_kwargs}
        **kwargs: {analysis_kwargs}

    """

    def __init__(
        self,
        data: AnnData,
        pval: float = 0.01,
        importance_cutoff: Number = 0.5,
        layer_key: Optional[str] = None,
        tree_kwargs: Optional[Dict] = None,
        **kwargs,
    ):
        try:
            from lightgbm import LGBMRegressor
        except ImportError:
            raise ImportError("lightgbm is not installed, please try `pip install lightgbm`.")
        super().__init__(data, display_name="NMD marker", **kwargs)

        if not self.neighbors_exists:
            raise NeighborsNotFoundError("Run `find_neighbors` first before continue.")

        tree_kwargs_ = {"n_jobs": -1, "random_state": 0, "importance_type": "gain"}
        if tree_kwargs is not None:
            for k, v in tree_kwargs.items():
                tree_kwargs_[k] = v

        markers = self.data.var[self.marker_key]
        neighbors = read_neighbors(self.data.obs, self.neighbors_key)
        cent_exp = read_exp(self.data, layer_key)
        # treat the neighbors as single cell
        # sum the expression
        neigh_exp = np.asarray([read_exp(self.data[n, :], layer_key).sum(1) for n in neighbors])
        results_data = []
        for i, y in enumerate(pbar_iter(cent_exp, desc="Neighbor-depedent markers",)):
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

        self.result = pd.DataFrame(
            data=results_data,
            columns=["marker", "neighbor_marker", "dependency", "corr", "pval"],
        )
