import pandas as pd
from anndata import AnnData
from spatialtis_core import CellCombs

from spatialtis.abc import AnalysisBase
from spatialtis.spatial.utils import NeighborsNotFoundError
from spatialtis.utils import doc
from spatialtis.utils.io import read_neighbors


@doc
class cell_interaction(AnalysisBase):
    """`Profiling cell-cell interaction <about/implementation.html#profiling-of-cell-cell-interaction>`_
    using permutation test

    Neighborhood analysis tells you the relationship between different type of cells

    - Association (1)
    - Avoidance (-1)
    - No relationship (0)

    This method is implemented in Rust, it executes in parallel automatically.

    Args:
        data: {adata}
        method: "pval" and "zscore" (Default: "pval")
        resample: Number of times to perform resample
        pval: {pval}
        order: If False, (Cell_A, Cell_B) and (Cell_B, Cell_A) are the same interaction (Default: False)
        **kwargs: {analysis_kwargs}

    .. seealso:: `spatial_enrichment_analysis <#spatialtis.plotting.spatial_enrichment_analysis>`_

    """

    def __init__(
            self,
            data: AnnData,
            method: str = "pval",
            resample: int = 1000,
            pval: float = 0.01,
            order: bool = False,
            **kwargs,
    ):
        if method == "pval":
            self.method = "pseudo p-value"
        else:
            self.method = "z-score"
        super().__init__(data, **kwargs)
        self.params = {"order": order}

        if not self.neighbors_exists:
            raise NeighborsNotFoundError("Run `find_neighbors` first before continue.")

        cc = CellCombs(self.cell_types, order)

        results_data = []
        for roi_name, roi_data in self.roi_iter(desc="Neighborhood analysis"):
            neighbors = read_neighbors(roi_data, self.neighbors_key)
            labels = roi_data[self.cell_id_key]
            cell_types = roi_data[self.cell_type_key]
            result = cc.bootstrap(cell_types, neighbors, labels, times=resample,
                                  pval=pval, method=method, ignore_self=True)
            for pairs in result:
                results_data.append([*roi_name, *pairs])

        df = pd.DataFrame(data=results_data, columns=self.exp_obs + ["type1", "type2", "value"])
        df = df.pivot_table(values="value", index=self.exp_obs, columns=["type1", "type2"])
        self.result = df
