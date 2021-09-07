from typing import Dict, Optional

import numpy as np
import pandas as pd
from anndata import AnnData
from lightgbm import LGBMRegressor
from scipy.stats import mannwhitneyu
from spatialtis_core import neighbor_components

from spatialtis.abc import AnalysisBase
from spatialtis.spatial.utils import NeighborsNotFoundError
from spatialtis.typing import Number
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
            importance_cutoff: Number = 0.5,
            pval: Number = 0.01,
            layer_key: Optional[str] = None,
            tree_kwargs: Optional[Dict] = None,
            **kwargs,
    ):
        super().__init__(data, **kwargs)

        if not self.neighbors_exists:
            raise NeighborsNotFoundError("Run `find_neighbors` first before continue.")

        if self.cell_type_key not in self.data.obs.keys():
            use_cell_type = False

        tree_kwargs_ = {"n_jobs": -1, "random_state": 0, "importance_type": "gain"}
        if tree_kwargs is not None:
            for k, v in tree_kwargs.items():
                tree_kwargs_[k] = v

        markers = self.data.var[self.marker_key]
        neighbors = read_neighbors(self.data, self.neighbors_key)
        labels = self.data.obs[self.cell_id_key]
        cell_types = self.data.obs[self.cell_type_key]
        col, comps = neighbor_components(neighbors, labels.tolist(), cell_types.tolist())
        neigh_comp = pd.DataFrame(data=comps,
                                  columns=col,
                                  index=pd.MultiIndex.from_frame(self.data.obs[[self.cell_type_key, self.cell_id_key]],
                                                                 names=['type', 'id'])
                                  )
        results_data = []
        # For markers in different cell types
        for t, x in neigh_comp.groupby(level=['type']):
            exp_ix = x.index.to_frame()['id']
            exp = read_exp(self.data[exp_ix, :], layer_key)
            for i, y in enumerate(exp):
                reg = LGBMRegressor(**tree_kwargs_).fit(x, y)
                weights = np.asarray(reg.feature_importances_)
                weights = weights / weights.sum()
                max_ix = np.argmax(weights)
                max_weight = weights[max_ix]
                max_type = col[max_ix]
                if max_weight > importance_cutoff:
                    nx = x.copy()
                    nx["exp"] = y  # add expression data to dataframe to allow cutting afterwards
                    at_neighbor = nx.iloc[:, max_ix] != 0  # cells with max_type at neighbors
                    at_neighbor_exp = nx[at_neighbor]["exp"].to_numpy()
                    non_at_neighbor_exp = nx[~at_neighbor]["exp"].to_numpy()
                    u, pvalue = mannwhitneyu(at_neighbor_exp, non_at_neighbor_exp)
                    if pvalue < pval:
                        log2_fc = np.log2(at_neighbor_exp.mean() / non_at_neighbor_exp.mean())
                        results_data.append([t, markers[i], max_type, max_weight, log2_fc, pvalue])
        self.result = pd.DataFrame(data=results_data,
                                   columns=["cell_type", "marker", "neighbor_type", "dependency", "log2_FC", "pval"],)
