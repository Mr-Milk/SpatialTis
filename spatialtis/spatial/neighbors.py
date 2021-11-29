from typing import Optional, List

import pandas as pd
from anndata import AnnData
from spatialtis_core import bbox_neighbors, multipoints_bbox, points_neighbors, spatial_weight

from spatialtis.abc import AnalysisBase
from spatialtis.typing import Number
from spatialtis.utils import col2adata_obs, doc, read_points, read_shapes, read_neighbors


@doc
class find_neighbors(AnalysisBase):
    """To `find the neighbors <../about/implementation.html#find-cell-neighbors>`_ of each cell

    KD-tree and Delaunay triangulation are used when cells are points

    R-tree is used when cells are polygons

    .. note::
        When :code:`method="kdtree"`, you can search neighbors within radius and/or by nearest-neighbors.
        If you specific :code:`r=30, k=5`, this will search within 30 while limited the number of neighbors
        to 5;

    Args:
        data: {adata}
        method: "kdtree", "rtree" and "delaunay", Default: "kdtree"
        r: The search radius
        k: The (minimum) number of nearest-neighbors
        scale: How much to scale each cell, only if cell has shape
        export_key: {export_key}
        **kwargs: {analysis_kwargs}

    Methods:
        spatial_weights: A generator that return spatial weight in CSR matrix (roi_name, spatial_weight_matrix)

    """

    def __init__(
        self,
        data: AnnData,
        method: Optional[str] = "kdtree",  # kdtree, delaunay, rtree,
        r: Optional[float] = None,
        k: Optional[int] = None,
        scale: Optional[Number] = None,
        export_key: Optional[str] = None,
        **kwargs,
    ):

        super().__init__(data,
                         method=method,
                         export_key=export_key,
                         **kwargs)

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
        export_key = self.neighbors_key if export_key is None else export_key
        col2adata_obs(neighbors_data, data, export_key)
        col2adata_obs(counts, data, f"{export_key}_count")
        self.stop_timer()

    @property
    def result(self):
        return self.data.obs[[self.cell_id_key, self.neighbors_key]]

    def spatial_weights(self):
        """A generator that return spatial weight in CSR matrix

        Returns:
            (roi_name, spatial_weight_matrix)

        """
        for roi_name, roi_data in self.roi_iter(disable_pbar=True):
            labels = roi_data[self.cell_id_key]
            neighbors = read_neighbors(roi_data, self.neighbors_key)
            yield roi_name, spatial_weight(neighbors, labels)
