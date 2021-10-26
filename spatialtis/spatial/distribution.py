from typing import Optional

import numpy as np
import pandas as pd
from anndata import AnnData
from spatialtis_core import points_bbox, spatial_distribution_pattern

from spatialtis.abc import AnalysisBase
from spatialtis.typing import Number, Tuple
from spatialtis.utils import doc, read_points


@doc
class cell_dispersion(AnalysisBase):
    """Cell distribution pattern

    There are three type of distribution pattern (0 if no cells)

    - Random (1)
    - Regular (2)
    - Cluster (3)

    Three methods are provided

    - `Index of Dispersion <../about/implementation.html#index-of-dispersion>`_
    - `Morisita’s index of dispersion <../about/implementation.html#morisitas-index-of-dispersion>`_
    - `Clark and Evans aggregation index <../about/implementation.html#clark-and-evans-aggregation-index>`_

    Notice that clark evans' index usually failed to detect local aggregation.

    +--------------------------------------+--------+---------+---------+
    |                                      | Random | Regular | Cluster |
    +======================================+========+=========+=========+
    | Index of dispersion: ID              | ID = 1 | ID < 1  | ID > 1  |
    +--------------------------------------+--------+---------+---------+
    | Morisita’s index of dispersion: I    | I = 1  |  I < 1  |  I > 1  |
    +--------------------------------------+--------+---------+---------+
    | Clark and Evans aggregation index: R | R = 1  |  R > 1  |  R < 1  |
    +--------------------------------------+--------+---------+---------+

    Args:
        data: {adata}
        method: "id", "morisita", and "clark_evans" (Default: "id")
        min_cells: The minimum number of the specific type of cells in a ROI to perform analysis
        pval: {pval}
        r: :code:`method="id"`, determine diameter of sample window, should be in [0, 1], default is 0.1
            this take 1/10 of the shortest side of the ROI as the diameter.
        resample: :code:`method="id"`, the number of random permutations to perform
        quad: :code:`method="morisita"`, {quad}
        rect_size: :code:`method="morisita"`, {rect_size}
        **kwargs: {analysis_kwargs}

    "quad" is quadratic statistic, it cuts a ROI into few rectangles, quad=(10,10) means the ROI will have 10*10 grid.

    """

    def __init__(
        self,
        data: AnnData,
        method: str = "id",
        min_cells: int = 10,
        pval: float = 0.01,
        r: Optional[Number] = None,
        resample: int = 1000,
        quad: Optional[Tuple[int, int]] = None,
        rect_size: Optional[Number] = None,
        **kwargs,
    ):
        display_method = {
            "id": "Index of dispersion",
            "morisita": "Morisita index",
            "clark_evans": "Clark evans index",
        }
        self.method = display_method[method]

        super().__init__(data, **kwargs)

        results_data = []
        for roi_name, roi_data in self.roi_iter(desc="Cell dispersion"):
            points = read_points(roi_data, self.centroid_key)
            bbox = points_bbox(points)
            new_df = pd.DataFrame(
                dict(points=points, cell_types=roi_data[self.cell_type_key])
            )
            points_collections = []
            cell_types = []
            for c, g in new_df.groupby("cell_types"):
                points_collections.append(g["points"])
                cell_types.append(c)
            result = spatial_distribution_pattern(
                points_collections,
                bbox,
                method=method,
                r=r,
                resample=resample,
                quad=quad,
                rect_side=rect_size,
                pval=pval,
                min_cells=min_cells,
            )
            for c, pattern in zip(cell_types, result):
                results_data.append([*roi_name, c, *pattern])
        results_data = pd.DataFrame(
            data=results_data,
            columns=self.exp_obs + ["cell_type", "index_value", "pval", "pattern"],
        )
        self.params = dict(exp_obs=self.exp_obs)
        self.result = results_data
