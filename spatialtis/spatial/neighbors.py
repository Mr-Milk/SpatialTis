from typing import Optional

import pandas as pd
from anndata import AnnData
from spatialtis_core import bbox_neighbors, multipoints_bbox, points_neighbors, points_bbox

from spatialtis.abc import AnalysisBase
from spatialtis.typing import Number
from spatialtis.utils import col2adata_obs, doc, read_points, read_shapes


@doc
class find_neighbors(AnalysisBase):
    """To `find the neighbors <../about/implementation.html#find-cell-neighbors>`_ of each cell

    KD-tree is used when cells are points

    R-tree is used when cells are polygons (use_shape=True)

    Args:
        data: {adata}
        expand: If the cell has shape, it means how much units to expand each cell;
                If the cell is point, it's the search radius;
        scale: How much to scale each cell, only if cell has shape
        count: Count the number of neighbors for each cell (Default: True)
        use_shape: Cell representation is shape instead of point (Default: False)
        **kwargs: {analysis_kwargs}

    """

    def __init__(
        self,
        data: AnnData,
        r: Optional[float] = None,
        k: Optional[int] = None,
        scale: Optional[Number] = None,
        method: Optional[str] = "kdtree",  # kdtree, delaunay, rtree,
        **kwargs,
    ):

        super().__init__(data, method=method, **kwargs)

        if method == "kdtree":
            if (r is None) & (k is None):
                raise ValueError(
                    "When search with KD-Tree, please specific"
                    "search radius `r` or number of neighbors `k`"
                )
        elif method == "rtree":
            if (r is None) & (scale is None):
                raise ValueError(
                    "When search with R-Tree, please specific"
                    "search radius `r` or scale ratio `scale`"
                )
        if scale is not None:
            if scale < 1:
                raise ValueError("Can't shrink cell, 'scale' must >= 1")

        if r is not None:
            if r < 0:
                raise ValueError("`r` must be greater than 0")

        if k is not None:
            if k < 0:
                raise ValueError("`k` must be greater than 0")
        # assign unique id to each cell, in case of someone cut the data afterwards
        # this ensure the analysis still work with non-integrated AnnData

        data.obs[self.cell_id_key] = [i for i in range(len(data.obs))]

        track_ix = []
        neighbors_data = []
        if method == "rtree":
            for roi_name, roi_data in self.roi_iter(desc="Find neighbors"):
                shapes = read_shapes(roi_data, self.shape_key)
                bbox = multipoints_bbox(shapes)
                labels = roi_data[self.cell_id_key]
                neighbors = bbox_neighbors(bbox, labels, expand=r, scale=scale)
                neighbors_data += neighbors
                track_ix += list(roi_data.index)
        else:
            for roi_name, roi_data in self.roi_iter(desc="Find neighbors"):
                points = read_points(roi_data, self.centroid_key)
                labels = roi_data[self.cell_id_key]
                neighbors = points_neighbors(points, labels, r=r, k=k, method=method)
                neighbors_data += neighbors
                track_ix += list(roi_data.index)
        counts = [len(i) for i in neighbors_data]
        neighbors_data = [str(i) for i in neighbors_data]
        neighbors_data = pd.Series(neighbors_data, index=track_ix)
        counts = pd.Series(counts, index=track_ix)
        col2adata_obs(neighbors_data, data, self.neighbors_key)
        col2adata_obs(counts, data, f"{self.neighbors_key}_count")
        self.stop_timer()

    @property
    def result(self):
        return self.data.obs[[self.cell_id_key, self.neighbors_key]]
