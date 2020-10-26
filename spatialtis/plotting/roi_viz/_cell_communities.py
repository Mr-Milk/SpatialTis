from ast import literal_eval
from collections import Counter
from typing import Dict, Optional

from anndata import AnnData

from spatialtis.config import CONFIG
from spatialtis.plotting.base import graph_plot, graph_plot_interactive
from spatialtis.utils import reuse_docstring


@reuse_docstring()
def cell_communities(
    adata: AnnData,
    query: Dict,
    min_cell: int = 10,
    use: str = "interactive",
    centroid_key: Optional[str] = None,
    community_key: Optional[str] = None,
    neighbors_key: Optional[str] = None,
    **kwargs,
):
    """(pyecharts, matplotlib) Visualize cell communities

    Args:
        adata: {adata_plotting}
        query: {query}
        min_cell: Only show communities with at least a number of cells
        use: Options are "interactive" and "_static". For big ROI, "interactive" is much faster using WebGL
        centroid_key: {centroid_key}
        community_key: {community_key}
        neighbors_key: {neighbors_key}
        **kwargs: Pass to `Graph plot (interactive) <plotting.html#spatialtis.plotting.graph_plot_interactive>`_
                  or `Graph plot <plotting.html#spatialtis.plotting.graph_plot>`_


    """
    if community_key is None:
        community_key = CONFIG.community_key
    if neighbors_key is None:
        neighbors_key = CONFIG.neighbors_key
    if centroid_key is None:
        centroid_key = CONFIG.CENTROID_KEY

    if neighbors_key not in adata.obs.keys():
        raise KeyError(
            "Neighbor key not found, export neighbors or specific your own key."
        )
    if community_key not in adata.obs.keys():
        raise KeyError(
            "Community key not found, run community or specific your own key."
        )

    df = adata.obs.query("&".join([f"({k}=='{v}')" for k, v in query.items()]))

    nodes_types = df[community_key].tolist()
    commus = []
    for commu, count in Counter(nodes_types).items():
        if count >= min_cell:
            commus.append(commu)

    df = df.reset_index(drop=True)
    xdf = df[df[community_key].isin(commus)]
    xdf = xdf.reset_index()

    nodes = [literal_eval(n) for n in xdf[centroid_key]]
    neighs = [literal_eval(n) for n in xdf[neighbors_key]]
    nodes_types = xdf[community_key]

    edges = []
    edges_types = []
    for i, n in zip(xdf.index, neighs):
        for x in n:
            new_x = xdf[xdf["index"] == x].index
            if len(new_x) == 1:
                new_x = new_x[0]
                if nodes_types[i] == nodes_types[new_x]:
                    edges.append((i, new_x))
                    edges_types.append(nodes_types[i])

    if use == "interactive":
        p = graph_plot_interactive(nodes, edges, edges_types=edges_types, **kwargs)
    else:
        p = graph_plot(nodes, edges, edges_types=edges_types, **kwargs)

    return p
