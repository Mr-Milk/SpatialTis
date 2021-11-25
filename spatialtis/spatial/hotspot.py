from typing import Optional, Tuple

import numpy as np
import pandas as pd
from anndata import AnnData
from scipy.spatial import cKDTree
from scipy.stats import norm
from spatialtis_core import getis_ord, points_bbox

from spatialtis.abc import AnalysisBase
from spatialtis.spatial.utils import QuadStats
from spatialtis.typing import Array
from spatialtis.utils import col2adata_obs, doc, read_points


@doc
class hotspot(AnalysisBase):
    """`Getis-ord hotspot detection <../about/implementation.html#hotspot-detection>`_

    Used to identify cells that cluster together.

    Args:
        data: {adata}
        selected_types: {selected_types}
        search_level: How deep the search level to reach
        quad: {quad}
        rect_size: {rect_size}
        pval: {pval}
        **kwargs: {analysis_kwargs}

    """

    def __init__(
        self,
        data: AnnData,
        selected_types: Optional[Array] = None,
        search_level: int = 3,
        quad: Optional[Tuple[int, int]] = None,
        rect_side: Optional[Tuple[float, float]] = None,
        pval: float = 0.01,
        **kwargs,
    ):
        super().__init__(data, **kwargs)
        self.export_key = "hotspot_all"
        if selected_types is not None:
            self.export_key = f"hotspot_{'_'.join(selected_types)}"
        else:
            selected_types = self.cell_types
        hotcells = []
        for roi_name, roi_data in self.roi_iter(desc="Hotspot analysis"):
            points = read_points(roi_data, self.centroid_key)
            bbox = points_bbox(points)
            for t, g in roi_data.groupby(self.cell_type_key):
                cells = read_points(g, self.centroid_key)
                if t in selected_types:
                    hots = getis_ord(
                        cells,
                        bbox,
                        search_level=search_level,
                        quad=quad,
                        rect_side=rect_side,
                        pval=pval,
                    )
                    hotcells.append(pd.Series(hots, index=g.index))

        result = pd.concat(hotcells)
        self.data.obs[self.export_key] = result
        # Cell map will leave blank if fill with None value
        self.data.obs[self.export_key].fillna("other", inplace=True)
        arr = self.data.obs[self.export_key].astype("category")
        arr = arr.cat.rename_categories({True: "hot", False: "cold", "other": "other"})
        self.data.obs[self.export_key] = arr
        # Call this to invoke the print
        col2adata_obs(self.data.obs[self.export_key], self.data, self.export_key)
        self.stop_timer()

    @property
    def result(self):
        return self.data.obs[self.exp_obs + [self.cell_type_key, self.export_key]]
