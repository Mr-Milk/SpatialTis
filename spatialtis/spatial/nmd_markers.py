import warnings
from typing import Dict, Optional

import numpy as np
import pandas as pd
from anndata import AnnData
from lightgbm import LGBMRegressor
from scipy.stats import spearmanr

from spatialtis.abc import AnalysisBase
from spatialtis.spatial.utils import NeighborsNotFoundError
from spatialtis.typing import Array, Number
from spatialtis.utils import doc
from spatialtis.utils.log import pbar_iter


@doc
class NMDMarkers(AnalysisBase):
    """Identify neighbor markers dependent marker

    Similar to `NCDMarkers`

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
            exp_std_cutoff: Number = 1.0,
            pval: float = 0.01,
            selected_markers: Optional[Array] = None,
            layers_key: Optional[str] = None,
            tree_kwargs: Optional[Dict] = None,
            **kwargs,
    ):
        super().__init__(data, task_name="NMDMarkers", **kwargs)

        if not self.neighbors_exists:
            raise NeighborsNotFoundError("Run `find_neighbors` first before continue.")

        tree_kwargs_ = {"n_jobs": -1, "random_state": 0, "importance_type": "gain"}
        if tree_kwargs is not None:
            for k, v in tree_kwargs.items():
                tree_kwargs_[k] = v

        results_data = []
        cent_cells, neigh_cells = self.get_neighbors_ix_pair()
        markers, (cent_exp, neigh_exp), _ = self.get_exp_matrix_fraction(
            markers=selected_markers,
            layers_key=layers_key,
            neighbors=(cent_cells, neigh_cells),
        )
        if len(markers) > 0:
            cent_exp = cent_exp.T

            for ix, m in enumerate(pbar_iter(markers, desc="NMD Markers", )):
                y = cent_exp[ix].copy()
                if np.std(y) > exp_std_cutoff:
                    reg = LGBMRegressor(**tree_kwargs_).fit(neigh_exp, y)
                    weights = np.asarray(reg.feature_importances_)
                    weights = weights / weights.sum()
                    max_ix = np.argmax(weights)
                    max_weight = weights[max_ix]
                    max_marker = markers[max_ix]
                    corr, pvalue = spearmanr(y, neigh_exp.T[max_ix])
                    if pvalue < pval:
                        results_data.append([m, max_marker, max_weight, corr, pvalue])
        else:
            warnings.warn("No markers left after filtering.")
        self.result = pd.DataFrame(
            data=results_data,
            columns=["marker", "neighbor_marker", "dependency", "corr", "pvalue"],
        )
