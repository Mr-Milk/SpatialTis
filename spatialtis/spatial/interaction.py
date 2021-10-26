import numpy as np
import pandas as pd
from anndata import AnnData
from spatialtis_core import CellCombs

from spatialtis.abc import AnalysisBase
from spatialtis.utils import NeighborsNotFoundError
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

    Args:
        data: {adata}
        method: "pval" and "zscore" (Default: "pval")
        resample: Number of times to perform resample
        pval: {pval}
        order: If False, (Cell_A, Cell_B) and (Cell_B, Cell_A) are the same interaction (Default: False)
        **kwargs: {analysis_kwargs}

    .. seealso:: :class:`spatialtis.spatial_enrichment`

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
        roi_tracker = []
        repeat_time = 0
        for roi_name, roi_data in self.roi_iter(desc="Cell interaction"):
            neighbors = read_neighbors(roi_data, self.neighbors_key)
            labels = roi_data[self.cell_id_key]
            cell_types = roi_data[self.cell_type_key]
            result = cc.bootstrap(
                cell_types,
                neighbors,
                labels,
                times=resample,
                pval=pval,
                method=method,
                ignore_self=True,
            )
            results_data += result
            roi_tracker += [roi_name]
            repeat_time = len(result)
        df = pd.DataFrame(data=results_data, columns=["type1", "type2", "value"])
        ix = pd.DataFrame(data=np.repeat(roi_tracker, repeat_time, axis=0), columns=self.exp_obs)
        df = pd.concat([ix, df], axis=1)
        df = df.pivot_table(values="value", index=self.exp_obs, columns=["type1", "type2"])

        self.result = df
