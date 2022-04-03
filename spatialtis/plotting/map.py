from itertools import cycle
from typing import List, Optional

import numpy as np
from anndata import AnnData
from matplotlib import pyplot as plt
from matplotlib.colors import to_hex
from milkviz import point_map, point_map3d, polygon_map
from natsort import natsorted
from scipy.sparse import issparse

from spatialtis import Config
from spatialtis.abc import AnalysisBase
from spatialtis.utils import doc, read_neighbors, read_points, read_shapes
from .utils import COLOR_POOL


@doc
def cell_map(
        data: AnnData,
        roi: str,
        use_shape: bool = False,
        selected_types: Optional[List] = None,
        masked_type_name: str = "Other",
        masked_type_color: str = "#d3d3d3",
        cell_type_key: Optional[str] = None,
        shape_key: Optional[str] = None,
        centroid_key: Optional[str] = None,
        roi_key: Optional[str] = None,
        **plot_options,
):
    """Visualize cells in ROI

    Args:
        data: {adata_plotting}
        roi: {roi}
        use_shape: Plot cell in polygon when shape data is available
        selected_types: {selected_types}
        masked_type_name: The name of the cell types not in selected_types
        masked_type_color: The color of the cell types not in selected_types
        cell_type_key: {cell_type_key}
        shape_key: {shape_key}
        centroid_key: {centroid_key}
        roi_key: {roi_key}
        **plot_options: Pass to `milkviz.point_map` or `milkviz.point_map3d` or `milkviz.polygon_map`

    """
    ab = AnalysisBase(data,
                      cell_type_key=cell_type_key,
                      shape_key=shape_key,
                      centroid_key=centroid_key,
                      roi_key=roi_key,
                      verbose=False)
    # cell_type_key = Config.cell_type_key if cell_type_key is None else cell_type_key
    # shape_key = Config.shape_key if shape_key is None else shape_key
    # centroid_key = Config.centroid_key if centroid_key is None else centroid_key
    # roi_key = Config.roi_key if roi_key is None else roi_key
    masked_type_color = to_hex(masked_type_color, keep_alpha=True)

    if ab.has_cell_type:
        all_cell_types = ab.cell_types
        color_mapper = dict(zip(all_cell_types, cycle(COLOR_POOL)))
        color_mapper[masked_type_name] = masked_type_color

    iter_data = data.obs.copy()
    points = ab.get_points()
    if len(points[0]) == 3:
        ab.dimension = 3
    iter_data['__spatial_centroid'] = points
    roi_info = iter_data[iter_data[ab.roi_key] == roi]

    if len(roi_info) == 0:
        raise ValueError(f"ROI not exist, roi = {roi}")
    cell_types = roi_info[ab.cell_type_key] if ab.has_cell_type else None

    internal_kwargs = dict(legend_title="Cell type")

    if selected_types is not None:
        utypes = np.unique(selected_types)
        cell_mask = cell_types.isin(utypes)
        cell_types = cell_types.to_numpy()
        cell_types[~cell_mask] = np.unique(masked_type_name)
        internal_kwargs["colors"] = [color_mapper.get(c) for c in cell_types]

    internal_kwargs = {**internal_kwargs, **plot_options}
    if use_shape:
        polygons = read_shapes(roi_info, ab.shape_key)
        return polygon_map(polygons, types=cell_types, **internal_kwargs)
    else:
        cells = np.array(roi_info['__spatial_centroid'].to_list())
        if ab.dimension == 2:
            x, y = cells[:, 0], cells[:, 1]
            return point_map(x, y, types=cell_types, **internal_kwargs)
        else:
            x, y, z = cells[:, 0], cells[:, 1], cells[:, 2]
            return point_map3d(x, y, z, types=cell_types, **internal_kwargs)


@doc
def expression_map(
        data: AnnData,
        roi: str,
        marker: str,
        use_shape: bool = False,
        selected_types: Optional[List] = None,
        cell_type_key: Optional[str] = None,
        marker_key: Optional[str] = None,
        shape_key: Optional[str] = None,
        centroid_key: Optional[str] = None,
        roi_key: Optional[str] = None,
        **plot_options,
):
    """Visualize marker expression in ROI

    Args:
        data: {adata_plotting}
        roi: {roi}
        marker:
        use_shape:
        marker_key:
        shape_key:
        centroid_key:
        roi_key:
        **plot_options:

    Returns:

    """
    ab = AnalysisBase(data,
                      # cell_type_key=cell_type_key,
                      shape_key=shape_key,
                      centroid_key=centroid_key,
                      roi_key=roi_key,
                      marker_key=marker_key,
                      verbose=False)
    # marker_key = Config.marker_key if marker_key is None else marker_key
    # shape_key = Config.shape_key if shape_key is None else shape_key
    # centroid_key = Config.centroid_key if centroid_key is None else centroid_key
    # roi_key = Config.roi_key if roi_key is None else roi_key

    internal_kwargs = dict(legend_title="expression")

    iter_data = data.obs.copy()
    points = ab.get_points()
    if len(points[0]) == 3:
        ab.dimension = 3
    iter_data['__spatial_centroid'] = points
    roi_info = iter_data[iter_data[ab.roi_key] == roi]

    roi_selector = data.obs[ab.roi_key] == roi

    if len(roi_info) == 0:
        raise ValueError(f"ROI not exist, roi = {roi}")
    if ab.marker_key is None:
        marker_v = data[roi_selector, data.var.index == marker].X.copy()
    else:
        marker_v = data[roi_selector, data.var[ab.marker_key] == marker].X.copy()
    if issparse(marker_v):
        marker_v = marker_v.A
    marker_v = marker_v.flatten()

    cell_mask = None
    if selected_types is not None:
        cell_type_key = Config.cell_type_key if cell_type_key is None else cell_type_key
        cell_types = roi_info[cell_type_key]
        utypes = np.unique(selected_types)
        cell_mask = cell_types.isin(utypes)

    internal_kwargs = {**internal_kwargs, **plot_options}
    if use_shape:
        polygons = read_shapes(roi_info, shape_key)
        if cell_mask is not None:
            polygons = np.asarray(polygons)[cell_mask]
            marker_v = marker_v[cell_mask]
        ax = polygon_map(polygons, values=marker_v, **internal_kwargs)
    else:
        cells = np.array(roi_info['__spatial_centroid'].to_list())
        if cell_mask is not None:
            cells = cells[cell_mask]
            marker_v = marker_v[cell_mask]
        if ab.dimension == 2:
            x, y = cells[:, 0], cells[:, 1]
            ax = point_map(x, y, values=marker_v, **internal_kwargs)
        else:
            x, y, z = cells[:, 0], cells[:, 1], cells[:, 2]
            ax = point_map3d(x, y, z, values=marker_v, **internal_kwargs)
    plt.title(f"{marker}")
    return ax


@doc
def neighbors_map(
        data: AnnData,
        roi: str,
        cell_type_key: Optional[str] = None,
        centroid_key: Optional[str] = None,
        roi_key: Optional[str] = None,
        **plot_options,
):
    """Visualize neighbors network built in a ROI

    Args:
        data: {adata_plotting}
        roi: {roi}
        cell_type_key: {cell_type_key}
        centroid_key: {centroid_key}
        roi_key: {roi_key}
        **plot_options:

    Returns:

    """
    # cell_type_key = Config.cell_type_key if cell_type_key is None else cell_type_key
    # centroid_key = Config.centroid_key if centroid_key is None else centroid_key
    # roi_key = Config.roi_key if roi_key is None else roi_key
    ab = AnalysisBase(data, cell_type_key=cell_type_key, centroid_key=centroid_key, roi_key=roi_key)

    iter_data = data.obs.copy()
    points = ab.get_points()
    if len(points[0]) == 3:
        raise NotImplementedError("Does not support 3D neighbor map")
    iter_data['__spatial_centroid'] = points
    roi_info = iter_data[iter_data[ab.roi_key] == roi]

    if len(roi_info) == 0:
        raise ValueError(f"ROI not exist, roi = {roi}")
    cell_types = roi_info[ab.cell_type_key] if ab.has_cell_type else None

    internal_kwargs = dict(legend_title="Cell type", **plot_options)

    cells = np.array(roi_info['__spatial_centroid'].to_list())
    x, y = cells[:, 0], cells[:, 1]
    neighbors = read_neighbors(roi_info, "cell_neighbors")
    labels = roi_info["cell_id"].astype(int)
    nmin = labels.min()
    links = []
    for l, neigh in zip(labels, neighbors):
        for n in neigh:
            if n > l:
                links.append((n - nmin, l - nmin))

    return point_map(x, y, types=cell_types, links=links, **internal_kwargs)
