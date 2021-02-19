from ast import literal_eval
from typing import Dict, Optional

from anndata import AnnData

from spatialtis import CONFIG
from spatialtis.plotting.base import expression_map_3d, expression_map_static
from spatialtis.utils import doc


@doc
def expression_map(
    data: AnnData,
    roi: Dict,
    marker: str,
    use: str = "static",
    marker_key: Optional[str] = None,
    centroid_key: Optional[str] = None,
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
        **plot_options: Pass to :class:`spatialtis.plotting.base.expression_map_static` or
            :class:`spatialtis.plotting.base.expression_map_3d`

    """
    if marker_key is None:
        marker_key = CONFIG.MARKER_KEY
    if centroid_key is None:
        centroid_key = CONFIG.CENTROID_KEY

    locations = data.obs.query(
        "&".join([f"({k}=='{v}')" for k, v in roi.items()])
    ).copy()[centroid_key]
    loc_ix = locations.index
    need_eval = isinstance(locations[0], str)
    if need_eval:
        locations = [literal_eval(loc) for loc in locations]

    marker_ix = data.var[marker_key].tolist().index(marker)
    expression = data[loc_ix].X.T[marker_ix]
    plot_options["saved_name"] = "expression_map_" + ",".join(
        [f"{k}={v}" for k, v in roi.items()]
    )
    if use == "interactive":
        return expression_map_3d(locations, expression, **plot_options)
    else:
        return expression_map_static(locations, expression, **plot_options)
