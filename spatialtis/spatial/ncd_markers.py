from collections import Counter
from typing import Dict, Optional

import numpy as np
import pandas as pd
from anndata import AnnData
from lightgbm import LGBMRegressor
from scipy.stats import spearmanr
from tqdm import tqdm

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
        use_cell_type: bool = False,
        corr_cutoff: Number = 0.5,
        importance_cutoff: Number = 0.6,
        exp_std_cutoff: Number = 1.0,
        pval: float = 0.01,
        selected_markers: Optional[Array] = None,
        layers_key: Optional[str] = None,
        tree_kwargs: Optional[Dict] = None,
        **kwargs,
    ):
        super().__init__(data, task_name="NCDMarkers", **kwargs)

        if not self.neighbors_exists:
            raise NeighborsNotFoundError("Run `find_neighbors` first before continue.")

        if self.cell_type_key not in self.data.obs.keys():
            use_cell_type = False

        tree_kwargs_ = {"n_jobs": -1, "random_state": 0, "importance_type": "gain"}
        if tree_kwargs is not None:
            for k, v in tree_kwargs.items():
                tree_kwargs_[k] = v

        if use_cell_type:
            result_data = []

            neighbors = self.get_neighbors_ix_map()
            cent_cells = list(neighbors.keys())
            cent_type = self.data.obs[self.cell_type_key][cent_cells]
            for t, g in tqdm(
                pd.DataFrame(
                    {"cent_cell": cent_cells, "cent_type": cent_type,}
                ).groupby("cent_type"),
                **CONFIG.pbar(desc="NCD Markers"),
            ):
                cents = g["cent_cell"].values
                markers, exp_matrix, _ = self.get_exp_matrix_fraction(
                    markers=selected_markers,
                    layers_key=layers_key,
                    neighbors_ix=cents,
                    std=exp_std_cutoff,
                )
                if len(markers) > 0:
                    neigh_types = pd.DataFrame(
                        [
                            Counter(self.data.obs[self.cell_type_key][neighbors[i]])
                            for i in cents
                        ]
                    )
                    cent_exp = exp_matrix.T
                    for ix, m in enumerate(markers):
                        y = cent_exp[ix].copy()
                        max_type, max_weight, corr, pvalue = max_contri_marker(
                            neigh_types, y, tree_kwargs_
                        )
                        if (
                            (pvalue < pval)
                            & (abs(corr) > corr_cutoff)
                            & (max_weight > importance_cutoff)
                        ):
                            result_data.append(
                                [t, m, max_type, max_weight, corr, pvalue]
                            )

            self.result = pd.DataFrame(
                data=result_data,
                columns=[
                    "cell_type",
                    "marker",
                    "neighbor_type",
                    "dependency",
                    "corr",
                    "pvalue",
                ],
            )

        else:
            neighbors = self.get_neighbors_ix_map()
            cent_cells = list(neighbors.keys())
            markers, exp_matrix, data = self.get_exp_matrix_fraction(
                markers=selected_markers,
                layers_key=layers_key,
                neighbors_ix=cent_cells,
                std=exp_std_cutoff,
            )
            if len(markers) > 0:
                neigh_types = pd.DataFrame(
                    [
                        Counter(self.data.obs[self.cell_type_key][neighbors[i]])
                        for i in cent_cells
                    ]
                )

                results_data = []
                cent_exp = exp_matrix.T
                for ix, m in enumerate(
                    tqdm(markers, **CONFIG.pbar(desc="NCD Markers"))
                ):
                    y = cent_exp[ix].copy()
                    max_type, max_weight, corr, pvalue = max_contri_marker(
                        neigh_types, y, tree_kwargs_
                    )
                    if (
                        (pvalue < pval)
                        & (abs(corr) > corr_cutoff)
                        & (max_weight > importance_cutoff)
                    ):
                        results_data.append([m, max_type, max_weight, corr, pvalue])

            self.result = pd.DataFrame(
                data=results_data,
                columns=["marker", "neighbor_type", "dependency", "corr", "pvalue"],
            )


def max_contri_marker(x, y, tree_kw):
    reg = LGBMRegressor(**tree_kw).fit(x, y)
    weights = np.asarray(reg.feature_importances_)
    weights = (weights - weights.min()) / (weights.max() - weights.min())
    max_ix = np.argmax(weights)
    max_weight = weights[max_ix]
    max_type = x.columns[max_ix]
    corr, pvalue = spearmanr(y, x[max_type])
    return max_type, max_weight, corr, pvalue
