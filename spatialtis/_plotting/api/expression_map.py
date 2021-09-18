from ast import literal_eval
from typing import Dict, Optional

import numpy as np
from anndata import AnnData
from milkviz import point_map
from spatialtis import Config
from spatialtis._plotting.api.utils import query_df
from spatialtis._plotting.base import expression_map_3d, expression_map_static
from spatialtis.utils import doc, read_points, read_exp


@doc
def expression_map(
    data: AnnData,
    roi: Dict,
    marker: str,
    marker_key: Optional[str] = None,
    centroid_key: Optional[str] = None,
    roi_key: Optional[str] = None,
    **plot_options,
):
    """Visualize cell expression in ROI

    Args:
        data: {adata_plotting}
        roi: {roi}
        marker: Which marker to visualize
        use: "static" or "interactive" (Default: "static")
        marker_key: {marker_key}
        centroid_key: {centroid_key}
        **plot_options: Pass to :class:`spatialtis._plotting.base.expression_map_static` or
            :class:`spatialtis._plotting.base.expression_map_3d`

    """
    marker_key = Config.marker_key if marker_key is None else marker_key
    centroid_key = Config.centroid_key if centroid_key is None else centroid_key
    roi_key = Config.roi_key if roi_key is None else roi_key

    roi_mask = data.obs[roi_key] == roi
    marker_mask = (data.var.index == marker) if marker_key is None else (data.var.index == marker)
    roi_info = data.obs[roi_mask]
    cells = np.array(read_points(roi_info, centroid_key))
    x, y = cells[:, 0], cells[:, 1]
    exp = data[roi_mask, marker_mask].X.T

    return point_map(x, y, values=exp, **plot_options)
