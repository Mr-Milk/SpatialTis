from collections import Counter
from typing import Dict, Optional

import numpy as np
import pandas as pd
from anndata import AnnData
from scipy.stats import spearmanr
from tqdm import tqdm
from xgboost import XGBRegressor

from spatialtis import CONFIG
from spatialtis.abc import AnalysisBase
from spatialtis.spatial.utils import NeighborsNotFoundError
from spatialtis.typing import Array, Number
from spatialtis.utils import doc


@doc
class NCDMarkers(AnalysisBase):
    """Identify neighbor cells dependent marker

    This method tells you the dependency and correlation between markers and its neighbor cell type.
    The dependency is calculated by building a gradiant boosting tree (in here XGBoost) to determine
    the feature importance. And the the spearman correlation is calculated.

    A reasonable std cutoff should be set, the marker expression need to have certain degree of variance.

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
        exp_std_cutoff: Number = 2.0,
        pval: float = 0.01,
        selected_markers: Optional[Array] = None,
        layers_key: Optional[str] = None,
        tree_kwargs: Optional[Dict] = None,
        **kwargs,
    ):
        super().__init__(data, task_name="NCDMarkers", **kwargs)

        if not self.neighbors_exists:
            raise NeighborsNotFoundError("Run `find_neighbors` first before continue.")

        tree_kwargs_ = {"n_jobs": -1, "random_state": 0}
        if tree_kwargs is not None:
            for k, v in tree_kwargs.items():
                tree_kwargs_[k] = v

        markers, exp_matrix, data = self.get_exp_matrix(selected_markers, layers_key)
        cent_cell, neigh_cells = self.get_neighbors_ix()
        cent_markers_exp = exp_matrix[cent_cell].T
        neigh_types = pd.DataFrame(
            [Counter(i) for i in data.obs[self.cell_type_key][neigh_cells]]
        )

        results_data = []
        for ix, m in enumerate(tqdm(markers, **CONFIG.pbar(desc="NCD Markers"))):
            y = cent_markers_exp[ix].copy()
            if np.std(y) >= exp_std_cutoff:
                reg = XGBRegressor(**tree_kwargs_).fit(neigh_types, y)
                weights = reg.feature_importances_
                max_ix = np.argmax(weights)
                max_weight = weights[max_ix]
                max_type = neigh_types.columns[max_ix]
                corr, pvalue = spearmanr(y, neigh_types[max_type])
                if pvalue < pval:
                    results_data.append([m, max_type, max_weight, corr, pvalue])

        self.result = pd.DataFrame(
            data=results_data,
            columns=["marker", "neighbor_type", "dependency", "corr", "pvalue"],
        )
