from ast import literal_eval

import numpy as np
import pandas as pd
from anndata import AnnData
from scipy.stats import norm

from spatialtis.abc import AnalysisBase
from spatialtis.spatial.utils import NeighborsNotFoundError
from spatialtis.utils import doc
from spatialtis.utils.log import pbar_iter


@doc
class neighborhood_analysis(AnalysisBase):
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
        super().__init__(data, task_name="neighborhood_analysis", **kwargs)
        self.params = {"order": order}

        try:
            import neighborhood_analysis as na
        except ImportError:
            raise ImportError(
                "Package not found, try `pip install neighborhood_analysis`."
            )

        if not self.neighbors_exists:
            raise NeighborsNotFoundError("Run `find_neighbors` first before continue.")

        need_eval = self.is_col_str(self.neighbors_key)
        cc = na.CellCombs(self.cell_types, order)

        results_data = []
        for name, g in pbar_iter(
                data.obs.groupby(self.exp_obs), desc="Neighborhood analysis"
        ):
            if isinstance(name, str):
                name = [name]
            if need_eval:
                neighbors = [literal_eval(n) for n in g[self.neighbors_key]]
            else:
                neighbors = [n for n in g[self.neighbors_key]]

            # here we need to move the neighbors from its real ix to starting point 0
            neighbors_min = np.asarray(
                [i for n in neighbors for i in n], dtype=int
            ).min()
            neighbors = [(np.asarray(n) - neighbors_min).tolist() for n in neighbors]

            cell_types = g[self.cell_type_key]
            result = cc.bootstrap(
                cell_types, neighbors, resample, pval, method, ignore_self=True
            )
            for k, v in result:
                results_data.append([*name, *k, v])

        df = pd.DataFrame(
            data=results_data, columns=self.exp_obs + ["type1", "type2", "value"]
        )

        if method == "pval":
            df.loc[:, ["value"]] = df["value"].astype(int)
        else:

            def sign(x):
                p = norm.sf(abs(x)) * 2
                if p < pval:
                    return np.sign(x)
                else:
                    return 0

            df.loc[:, ["value"]] = df["value"].apply(sign).astype(int)
        self.result = df
