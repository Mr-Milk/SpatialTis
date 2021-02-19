from itertools import combinations_with_replacement
from typing import Optional

import pandas as pd
from anndata import AnnData
from scipy.stats import pearsonr, spearmanr
from tqdm import tqdm

from spatialtis import CONFIG
from spatialtis.abc import AnalysisBase
from spatialtis.spatial.utils import NeighborsNotFoundError
from spatialtis.typing import Array, Number
from spatialtis.utils import doc


@doc
class spatial_co_expression(AnalysisBase):
    """Identifying spatial co-expression markers using correlation

    The correlation is calculated within pairs of neighbor cells

    Args:
        data: {adata}
        method: "spearman" or "pearson" (Default: "spearman")
        selected_markers: {selected_markers}
        layers_key: {layers_key}
        pval: {pval}

    """

    def __init__(
        self,
        data: AnnData,
        method: str = "spearman",
        selected_markers: Optional[Array] = None,
        layers_key: Optional[str] = None,
        pval: Number = 0.01,
        **kwargs,
    ):
        if method == "spearman":
            self.method = "spearman correlation"
            corr_func = spearmanr
        elif method == "pearson":
            self.method = "pearson correlation"
            corr_func = pearsonr
        else:
            raise ValueError("Available options are `spearman` and `pearson`.")
        super().__init__(data, task_name="spatial_co_expression", **kwargs)

        if not self.neighbors_exists:
            raise NeighborsNotFoundError("Run `find_neighbors` first before continue.")

        markers, exp_matrix, _ = self.get_exp_matrix(selected_markers, layers_key)
        cent_cells, neigh_cells = self.get_neighbors_ix()
        cent_exp = exp_matrix[cent_cells].T
        neigh_exp = exp_matrix[neigh_cells].T
        marker_ix = [i for i in combinations_with_replacement(range(len(markers)), 2)]

        result_data = []
        for p1, p2 in tqdm(marker_ix, **CONFIG.pbar(desc="co-expression")):
            corr_value, pvalue = corr_func(cent_exp[p1], neigh_exp[p2])
            if pvalue < pval:
                result_data.append([markers[p1], markers[p2], corr_value, pvalue])

        self.result = pd.DataFrame(
            result_data, columns=["marker1", "marker2", "corr", "pvalue"]
        )
