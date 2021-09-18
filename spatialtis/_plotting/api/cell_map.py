import warnings
from ast import literal_eval
from typing import Optional

import numpy as np
from anndata import AnnData
from milkviz import point_map, polygon_map

from spatialtis import Config
from spatialtis._plotting.base.cell_map import cell_map_interactive, cell_map_static
from spatialtis.typing import Array
from spatialtis.utils import doc, read_shapes, read_points


@doc
def cell_map(
    data: AnnData,
    roi: str,
    use_shape: bool = False,
    selected_types: Optional[Array] = None,
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
        use: "static" or "interactive" (Default: "static")
        selected_types: {selected_types}
        cell_type_key: {cell_type_key}
        shape_key: {shape_key}
        centroid_key: {centroid_key}
        roi_key: {roi_key}
        **plot_options: Pass to :class:`spatialtis._plotting.base.cell_map_static` or
            :class:`spatialtis._plotting.base.cell_map_interactive`

    """
    cell_type_key = Config.cell_type_key if cell_type_key is None else cell_type_key
    shape_key = Config.SHAPE_KEY if shape_key is None else shape_key
    centroid_key = Config.centroid_key if centroid_key is None else centroid_key
    roi_key = Config.roi_key if roi_key is None else roi_key

    roi_info = data.obs[data.obs[roi_key] == roi]
    cell_types = roi_info[cell_type_key]

    if use_shape:
        polygons = read_shapes(roi_info, shape_key)
        return polygon_map(polygons, cell_types, **plot_options)
    else:
        cells = np.array(read_points(roi_info, centroid_key))
        x, y = cells[:, 0], cells[:, 1]
        return point_map(x, y, types=cell_types, **plot_options)
