from itertools import combinations_with_replacement
from typing import Optional

import numpy as np
import pandas as pd
from anndata import AnnData
from spatialtis_core import fast_corr

from spatialtis.abc import AnalysisBase, neighbors_pairs
from spatialtis.typing import Number
from spatialtis.utils import doc, pbar_iter, NeighborsNotFoundError
from spatialtis.utils.io import read_exp, read_neighbors

DESCRIPTION = "co-expression"


@doc
class spatial_coexp(AnalysisBase):
    """Identifying spatial co-expression markers using correlation

    The correlation is calculated within pairs of neighbor cells

    Args:
        data: {adata}
        method: "spearman" or "pearson" (Default: "spearman")
        use_cel_type: Whether to use cell type information
        selected_markers: {selected_markers}
        layer_key: {layer_key}
        corr_thresh: The minimum correlation value to store the result,
        **kwargs: {analysis_kwargs}

    """

    def __init__(
            self,
            data: AnnData,
            method: str = "spearman",
            use_cell_type: bool = False,
            layer_key: Optional[str] = None,
            corr_thresh: Number = 0.5,
            **kwargs,
    ):
        if method == "spearman":
            self.method = "spearman correlation"
        elif method == "pearson":
            self.method = "pearson correlation"
        else:
            raise ValueError("Available options are `spearman` and `pearson`.")
        super().__init__(data, display_name="Spatial Co-expression", **kwargs)
        if not self.neighbors_exists:
            raise NeighborsNotFoundError("Run `find_neighbors` first before continue.")

        if self.cell_type_key not in self.data.obs.keys():
            use_cell_type = False

        self.params = {"use_cell_type": use_cell_type}

        neighbors = read_neighbors(self.data.obs, self.neighbors_key)
        labels = self.data.obs[self.cell_id_key]
        pairs = neighbors_pairs(labels, neighbors)
        markers = self.markers_col
        markers_combs = [(x, y) for x, y in combinations_with_replacement(markers, 2)]

        if use_cell_type:
            pairs_pool = {}
            pairs_order = {}
            type_pairs = self.data.obs[self.cell_type_key][pairs.ravel()]\
                .to_numpy().reshape(
                pairs.shape
            )
            for ix in range(type_pairs.shape[1]):
                i = type_pairs[:, ix]
                c = frozenset(i)
                if pairs_pool.get(c, 0) == 0:
                    pairs_pool[c] = 1
                    pairs_order[c] = i
                else:
                    if (pairs_order[c] != i).any():
                        pairs[:, ix] = pairs[:, ix][::-1]
                        type_pairs[:, ix] = type_pairs[:, ix][::-1]
            types = pd.DataFrame(
                np.vstack([pairs, type_pairs]).T, columns=["p1", "p2", "c1", "c2"]
            )

            data_collector = []
            for (t1, t2), df in pbar_iter(types.groupby(["c1", "c2"]), desc=DESCRIPTION):
                exp1 = read_exp(self.data[df["p1"].to_numpy(dtype=int), :], dtype=np.float, layer_key=layer_key)
                exp2 = read_exp(self.data[df["p2"].to_numpy(dtype=int), :], dtype=np.float, layer_key=layer_key)

                r = fast_corr(exp1, exp2, method=method)
                d = pd.DataFrame(markers_combs, columns=['marker1', 'marker2'])
                d['cell1'] = t1
                d['cell2'] = t2
                d['corr'] = r
                d = d[(d['corr'] > corr_thresh) | (d['corr'] < -corr_thresh)]
                data_collector.append(d)
            d = pd.concat(data_collector)
            self.result = d.sort_values('corr', ascending=False) \
                .reset_index(drop=True)

        else:
            exp1 = read_exp(self.data[pairs[0], :], dtype=np.float, layer_key=layer_key)
            exp2 = read_exp(self.data[pairs[1], :], dtype=np.float, layer_key=layer_key)
            r = fast_corr(exp1, exp2, method=method)
            d = pd.DataFrame(markers_combs, columns=['marker1', 'marker2'])
            d['corr'] = r
            self.result = d[(d['corr'] > corr_thresh) | (d['corr'] < -corr_thresh)] \
                .sort_values('corr', ascending=False) \
                .reset_index(drop=True)
