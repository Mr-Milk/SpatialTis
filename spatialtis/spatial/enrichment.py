import warnings
from ast import literal_eval
from typing import Optional

import numpy as np
import pandas as pd
from anndata import AnnData
from scipy.stats import norm
from tqdm import tqdm

from spatialtis.abc import AnalysisBase
from spatialtis.config import CONFIG
from spatialtis.spatial.utils import NeighborsNotFoundError
from spatialtis.typing import Array
from spatialtis.utils import doc


@doc
class spatial_enrichment_analysis(AnalysisBase):
    """`Profiling markers spatial enrichment <about/implementation.html#profiling-of-markers-co-expression>`_
    using permutation test

    Similar to neighborhood analysis which tells you the relationship between different type of cells.
    This analysis tells you the spatial relationship between markers.

    This method is implemented in Rust, it executes in parallel automatically.

    Args:
        data: {adata}
        threshold: The expression level to determine whether a marker is positive
        layers_key: {layers_key}
        selected_markers: {selected_markers}
        resample: Number of times to perform resample
        pval: {pval}
        order: If False, (Cell_A, Cell_B) and (Cell_B, Cell_A) are the same interaction (Default: False)
        **kwargs: {analysis_kwargs}

    .. seealso:: `neighborhood_analysis <#spatialtis.plotting.neighborhood_analysis>`_

    """

    def __init__(
        self,
        data: AnnData,
        threshold: Optional[float] = None,
        layers_key: Optional[str] = None,
        selected_markers: Optional[Array] = None,
        resample: int = 500,
        pval: float = 0.01,
        order: bool = False,
        **kwargs,
    ):
        super().__init__(data, task_name="spatial_enrichment_analysis", **kwargs)
        self.params = {"order": order}
        try:
            import neighborhood_analysis as na
        except ImportError:
            raise ImportError(
                "Package not found, try `pip install neighborhood_analysis`."
            )

        if not self.neighbors_exists:
            raise NeighborsNotFoundError("Run `find_neighbors` first before continue.")

        if (threshold is not None) & (layers_key is None):
            layers_key = f"gt_{threshold}"
            data.layers[layers_key] = (data.X.copy() >= threshold).astype(int)
        elif (threshold is None) & (layers_key is not None):
            warnings.warn(
                "You specific both threshold and layers_key, using user defined layers_key"
            )
        else:
            layers_key = f"mean_cut"
            data.layers[layers_key] = (data.X.copy() >= data.X.mean(axis=0)).astype(int)

        markers, _, data = self.get_exp_matrix(selected_markers, layers_key)
        need_eval = self.is_col_str(self.neighbors_key)
        results_data = []

        for name, roi in tqdm(
            data.obs.groupby(self.exp_obs),
            **CONFIG.pbar(desc="Spatial enrichment analysis"),
        ):
            if isinstance(name, str):
                name = [name]
            if need_eval:
                neighbors = [literal_eval(n) for n in roi[self.neighbors_key]]
            else:
                neighbors = [n for n in roi[self.neighbors_key]]
            matrix = data[roi.index].layers[layers_key]

            for ix, x in enumerate(markers):
                x_status = [bool(i) for i in matrix[:, ix]]
                for iy, y in enumerate(markers):
                    if (not order) & (iy < ix):
                        pass
                    else:
                        y_status = [bool(i) for i in matrix[:, iy]]
                        z = na.comb_bootstrap(
                            x_status,
                            y_status,
                            neighbors,
                            times=resample,
                            ignore_self=False,
                        )
                        results_data.append([*name, x, y, z])
        df = pd.DataFrame(
            data=results_data, columns=self.exp_obs + ["marker1", "marker2", "value"]
        )

        def sign(x):
            p = norm.sf(abs(x))
            if p < pval:
                return np.sign(x)
            else:
                return 0

        df.loc[:, ["value"]] = df["value"].apply(sign).astype(int)
        self.result = df
