import warnings
from ast import literal_eval
from typing import Optional

from anndata import AnnData

from spatialtis import Config
from spatialtis.plotting.base.cell_map import cell_map_interactive, cell_map_static
from spatialtis.typing import Array
from spatialtis.utils import doc


@doc
def cell_map(
        data: AnnData,
        roi: str,
        use_shape: bool = False,
        selected_types: Optional[Array] = None,
        use: str = "static",
        cell_type_key: Optional[str] = None,
        shape_key: Optional[str] = None,
        centroid_key: Optional[str] = None,
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
        **plot_options: Pass to :class:`spatialtis.plotting.base.cell_map_static` or
            :class:`spatialtis.plotting.base.cell_map_interactive`

    """
    if cell_type_key is None:
        cell_type_key = Config.cell_type_key
    if shape_key is None:
        shape_key = Config.SHAPE_KEY
    if centroid_key is None:
        centroid_key = Config.centroid_key

    if use_shape:
        if shape_key not in data.obs.keys():
            warnings.warn("Shape key not exist, try to resolve cell as point")
            data_key = centroid_key
        else:
            data_key = shape_key
    else:
        if centroid_key not in data.obs.keys():
            raise KeyError("Centroid key not exist")
        data_key = centroid_key

    data = data.obs[data.obs[Config.roi_key] == roi]
    need_eval = isinstance(data[data_key][0], str)
    plot_data = {}
    for n, df in data.groupby(cell_type_key):
        if need_eval:
            cells = [literal_eval(c) for c in df[data_key]]
        else:
            cells = [c for c in df[data_key]]
        plot_data[n] = cells
    plot_options["saved_name"] = "cell_map_" + roi
    if data_key == centroid_key:
        if use == "interactive":
            return cell_map_interactive(
                points=plot_data, selected_types=selected_types, **plot_options
            )
        else:
            return cell_map_static(
                points=plot_data, selected_types=selected_types, **plot_options
            )
    else:
        if use == "interactive":
            return cell_map_interactive(
                shapes=plot_data, selected_types=selected_types, **plot_options
            )
        else:
            return cell_map_static(
                shapes=plot_data, selected_types=selected_types, **plot_options
            )
