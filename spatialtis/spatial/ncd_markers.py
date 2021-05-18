from typing import Dict, Optional

import numpy as np
import pandas as pd
from anndata import AnnData
from lightgbm import LGBMRegressor
from scipy.stats import mannwhitneyu

from spatialtis.abc import AnalysisBase
from spatialtis.spatial.utils import NeighborsNotFoundError, normalize
from spatialtis.typing import Array, Number
from spatialtis.utils import doc
from spatialtis.utils.log import pbar_iter

try:
    import neighborhood_analysis as na
except ImportError:
    raise ImportError("Try pip install neighborhood_analysis")


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
            importance_cutoff: Number = 0.5,
            exp_std_cutoff: Number = 1.0,
            pval: Number = 0.01,
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
            df = self.data.obs.set_index(self.neighbors_ix_key)
            cell_types = df[self.cell_type_key]
            cent_type = cell_types[cent_cells]

            ix, col, data = na.neighbor_components(
                neighbors, dict(zip(df.index, cell_types))
            )
            neigh_comp = pd.DataFrame(data=data, index=ix, columns=col)

            markers, exp_matrix, cut_data = self.get_exp_matrix_fraction(
                markers=selected_markers, layers_key=layers_key,
            )

            if len(markers) > 0:
                for t, g in pbar_iter(
                        pd.DataFrame(
                            {"cent_cell": cent_cells, "cent_type": cent_type, }
                        ).groupby("cent_type"),
                        desc="NCD Markers",
                ):
                    cents = g["cent_cell"].values
                    meta = (
                        cut_data.obs.reset_index(drop=True)
                            .reset_index()
                            .set_index(self.neighbors_ix_key)
                    )
                    exp_ix = meta.loc[cents]["index"].values
                    exp = exp_matrix[exp_ix]

                    neigh_types = neigh_comp.loc[cents]
                    cent_exp = exp.T
                    for ix, m in enumerate(markers):
                        y = cent_exp[ix].copy()
                        max_type, max_weight, lo2_fc, pvalue = max_contri_marker(
                            neigh_types,
                            y,
                            exp_std_cutoff,
                            tree_kwargs_,
                            importance_cutoff,
                            pval,
                        )
                        if (max_type is not None) & (max_type != t):
                            result_data.append(
                                [t, m, max_type, max_weight, lo2_fc, pvalue]
                            )

                self.result = pd.DataFrame(
                    data=result_data,
                    columns=[
                        "cell_type",
                        "marker",
                        "neighbor_type",
                        "dependency",
                        "log2_FC",
                        "pvalue",
                    ],
                )

        else:
            results_data = []
            neighbors = self.get_neighbors_ix_map()
            cent_cells = list(neighbors.keys())
            df = self.data.obs.set_index(self.neighbors_ix_key)
            markers, exp_matrix, data = self.get_exp_matrix_fraction(
                markers=selected_markers,
                layers_key=layers_key,
                neighbors_ix=cent_cells,
            )
            if len(markers) > 0:
                ix, col, data = na.neighbor_components(
                    neighbors, dict(zip(df.index, df[self.cell_type_key]))
                )
                neigh_comp = pd.DataFrame(data=data, index=ix, columns=col)
                neigh_types = neigh_comp.loc[cent_cells]

                cent_exp = exp_matrix.T
                for ix, m in enumerate(
                        pbar_iter(markers, desc="NCD Markers", )
                ):
                    y = cent_exp[ix].copy()
                    max_type, max_weight, log2_fc, pvalue = max_contri_marker(
                        neigh_types,
                        y,
                        exp_std_cutoff,
                        tree_kwargs_,
                        importance_cutoff,
                        pval,
                    )
                    if max_type is not None:
                        results_data.append([m, max_type, max_weight, log2_fc, pvalue])

            self.result = pd.DataFrame(
                data=results_data,
                columns=["marker", "neighbor_type", "dependency", "log2_FC", "pvalue"],
            )


def max_contri_marker(x, y, exp_std_cutoff, tree_kw, importance_cutoff, pval):
    if np.std(y) > exp_std_cutoff:
        reg = LGBMRegressor(**tree_kw).fit(x, y)
        weights = np.asarray(reg.feature_importances_)
        weights = weights / weights.sum()
        if weights.mean() < 1:
            max_ix = np.argmax(weights)
            max_weight = weights[max_ix]
            max_type = x.columns[max_ix]
            if max_weight > importance_cutoff:
                x = x.copy()
                x["y"] = y
                at_neighbor = x.iloc[:, max_ix] != 0
                at_neighbor_exp = normalize(x[at_neighbor]["y"].to_numpy())
                non_at_neighbor_exp = normalize(x[~at_neighbor]["y"].to_numpy())
                u, pvalue = mannwhitneyu(at_neighbor_exp, non_at_neighbor_exp)
                log2_fc = np.log2(at_neighbor_exp.mean() / non_at_neighbor_exp.mean())
                if pvalue < pval:
                    return max_type, max_weight, log2_fc, pvalue

    return None, None, None, None
