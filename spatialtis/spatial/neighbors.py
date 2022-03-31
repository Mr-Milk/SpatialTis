from typing import Optional

import pandas as pd
from anndata import AnnData
from spatialtis_core import bbox_neighbors, multipoints_bbox, points_neighbors, spatial_weight
from spatialtis_core.neighbors import bbox_neighbors_parallel, points_neighbors_parallel

from spatialtis.abc import AnalysisBase
from spatialtis.typing import Number
from spatialtis.utils import col2adata_obs, doc, read_shapes, read_neighbors


@doc
def find_neighbors(data: AnnData,
                   method: Optional[str] = "kdtree",  # kdtree, delaunay, rtree,
                   r: Optional[float] = None,
                   k: Optional[int] = None,
                   scale: Optional[Number] = None,
                   export_key: Optional[str] = None,
                   **kwargs, ):
    """To `find the neighbors <../about/implementation.html#find-cell-neighbors>`_ of each cell

    KD-tree and Delaunay triangulation are used when cells are points, by default search for k=5

    R-tree is used when cells are polygons, by default search for scale=1.4

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

    ab = AnalysisBase(data, method=method,
                      display_name="Find neighbors",
                      export_key=export_key, **kwargs)

    if method == "kdtree":
        if (r is None) & (k is None):
            k = 5
        if k is not None:
            if k < 0:
                raise ValueError("`k` must be greater than 0")

    elif method == "rtree":
        if (r is None) & (scale is None):
            scale = 1.4
        if scale is not None:
            if scale < 1:
                raise ValueError("Can't shrink cell, 'scale' must >= 1")

    if r is not None:
        if r < 0:
            raise ValueError("`r` must be greater than 0")

    # assign unique id to each cell, in case of someone cut the data afterwards
    # this ensure the analysis still work with non-integrated AnnData

    data.obs[ab.cell_id_key] = [i for i in range(len(data.obs))]

    track_ix = []
    neighbors_data = []
    if method == "rtree":
        if ab.dimension == 3:
            raise NotImplementedError("Don't support RTree with 3D data")

        bbox_collections = []
        labels_collections = []
        for roi_name, roi_data in ab.roi_iter(desc="Find neighbors"):
            shapes = read_shapes(roi_data, ab.shape_key)
            bbox = multipoints_bbox(shapes)
            labels = roi_data[ab.cell_id_key].values.tolist()
            bbox_collections.append(bbox)
            labels_collections.append(labels)
            # neighbors = bbox_neighbors(bbox, labels, expand=r, scale=scale)
            # neighbors_data += neighbors
            track_ix += list(roi_data.index)
        neighbors = bbox_neighbors_parallel(bbox_collections, labels_collections, expand=r, scale=scale)
        for n in neighbors:
            neighbors_data += n
    else:

        points_collections = []
        labels_collections = []
        for roi_name, roi_data, points in ab.roi_iter_with_points(desc="Find neighbors"):
            labels = roi_data[ab.cell_id_key].values.tolist()
            points_collections.append(points)
            labels_collections.append(labels)
            track_ix += list(roi_data.index)

        neighbors = points_neighbors_parallel(points_collections, labels_collections, r=r, k=k, method=method)
        for n in neighbors:
            neighbors_data += n

    counts = [len(i) for i in neighbors_data]
    neighbors_data = [str(i) for i in neighbors_data]
    neighbors_data = pd.Series(neighbors_data, index=track_ix)
    counts = pd.Series(counts, index=track_ix)
    export_key = ab.neighbors_key if export_key is None else export_key
    col2adata_obs(neighbors_data, data, export_key)
    col2adata_obs(counts, data, f"{export_key}_count")
    ab.stop_timer()


@doc
def spatial_weights(
        data: AnnData,
        **kwargs
):
    """A generator that return spatial weight in CSR matrix

    Returns:
        (roi_name, spatial_weight_matrix)

    Examples:
        >>> import spatialtis as st
        >>> weights_matrix_every_roi = [matrix for (roi, matrix) in st.spatial_weights(data)]

    """
    ab = AnalysisBase(data, **kwargs)
    ab.check_neighbors()
    for roi_name, roi_data in ab.roi_iter(disable_pbar=True):
        labels = roi_data[ab.cell_id_key]
        neighbors = read_neighbors(roi_data, ab.neighbors_key)
        yield roi_name, spatial_weight(neighbors, labels)
