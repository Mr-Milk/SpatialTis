from ast import literal_eval
from typing import Dict, Optional

from anndata import AnnData

from spatialtis.config import CONFIG
from spatialtis.plotting.base import graph_position_interactive, graph_position_static
from spatialtis.utils import doc


@doc
def neighbors_map(
    data: AnnData,
    roi: Dict,
    use: str = "static",
    cell_type_key: Optional[str] = None,
    centroid_key: Optional[str] = None,
    neighbors_key: Optional[str] = None,
    **plot_options,
):
    """Visualize cell neighbors in ROI

    Args:
        data: {adata_plotting}
        roi: {roi}
        use: "static" or "interactive" (Default: "static")
        cell_type_key: {cell_type_key}
        centroid_key: {centroid_key}
        neighbors_key: {neighbors_key}
        **plot_options: Pass to :class:`spatialtis.plotting.base.graph_position_static` or
            :class:`spatialtis.plotting.base.graph_position_interactive`

    {pyecharts_tips}

    """
    if cell_type_key is None:
        cell_type_key = CONFIG.CELL_TYPE_KEY
    if centroid_key is None:
        centroid_key = CONFIG.CENTROID_KEY
    if neighbors_key is None:
        neighbors_key = CONFIG.NEIGHBORS_KEY

    df = data.obs.query("&".join([f"({k}=='{v}')" for k, v in roi.items()])).copy()
    df = df.reset_index()

    need_eval_nodes = isinstance(df[centroid_key][0], str)
    if need_eval_nodes:
        nodes = [literal_eval(n) for n in df[centroid_key]]
    else:
        nodes = [n for n in df[centroid_key]]
    nodes_types = [str(n) for n in df[cell_type_key]]

    need_eval_neighs = isinstance(df[neighbors_key][0], str)
    if need_eval_neighs:
        neighs = [literal_eval(n) for n in df[neighbors_key]]
    else:
        neighs = [n for n in df[neighbors_key]]

    edges = []
    for i, n in zip(df.index, neighs):
        for x in n:
            edges.append((i, x))

    plot_options["saved_name"] = "neighbors_map_" + ",".join(
        [f"{k}={v}" for k, v in roi.items()]
    )
    if use == "interactive":
        return graph_position_interactive(
            nodes, edges, nodes_types=nodes_types, node_size=4, **plot_options
        )
    else:
        return graph_position_static(
            nodes, edges, nodes_types=nodes_types, node_size=3, **plot_options
        )
