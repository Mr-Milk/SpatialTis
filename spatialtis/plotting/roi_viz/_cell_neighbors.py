from ast import literal_eval
from typing import Dict, Optional

from anndata import AnnData

from spatialtis.config import CONFIG
from spatialtis.plotting.base import graph_plot, graph_plot_interactive
from spatialtis.utils import reuse_docstring


@reuse_docstring()
def cell_neighbors(
    adata: AnnData,
    query: Dict,
    use: str = "interactive",
    type_key: Optional[str] = None,
    centroid_key: Optional[str] = None,
    neighbors_key: Optional[str] = None,
    **kwargs,
):
    """(pyecharts) Visualize cell neighbors in ROI

    Args:
        adata: {adata_plotting}
        query: {query}
        use: Options are "interactive" and "_static"
        type_key: {type_key}
        centroid_key: {centroid_key}
        neighbors_key: {neighbors_key}
        **kwargs: Pass to `Graph plot (interactive) <plotting.html#spatialtis.plotting.graph_plot_interactive>`_
            or `Graph plot <plotting.html#spatialtis.plotting.graph_plot>`_

    """
    if type_key is None:
        type_key = CONFIG.CELL_TYPE_KEY
    if neighbors_key is None:
        neighbors_key = CONFIG.neighbors_key
    if centroid_key is None:
        centroid_key = CONFIG.CENTROID_KEY

    df = adata.obs.query("&".join([f"({k}=='{v}')" for k, v in query.items()]))
    df = df.reset_index()

    nodes = [literal_eval(n) for n in df[centroid_key]]
    nodes_types = [str(n) for n in df[type_key]]

    neighs = [literal_eval(n) for n in df[neighbors_key]]
    edges = []
    for i, n in zip(df.index, neighs):
        for x in n:
            edges.append((i, x))

    if use == "interactive":
        p = graph_plot_interactive(nodes, edges, nodes_types=nodes_types, node_size=3)
    else:
        p = graph_plot(nodes, edges, nodes_types=nodes_types, node_size=4)
    return p
