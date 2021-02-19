from ast import literal_eval
from collections import Counter
from typing import Dict, Optional

from anndata import AnnData

from spatialtis.config import ANALYSIS, CONFIG

from ...utils import doc
from ..base import graph_position_interactive, graph_position_static


@doc
def community_map(
    data: AnnData,
    roi: Dict,
    min_cells: int = 10,
    use: str = "static",
    community_key: Optional[str] = None,
    centroid_key: Optional[str] = None,
    neighbors_key: Optional[str] = None,
    **plot_options,
):
    """Visualize cell communities in ROI

    Args:
        data: {adata_plotting}
        roi: {roi}
        min_cells: Show communities contain more than a number of cells
        use: "static" or "interactive" (Default: "static")
        community_key: {community_key}
        centroid_key: {centroid_key}
        neighbors_key: {neighbors_key}
        **plot_options: Pass to :class:`spatialtis.plotting.base.graph_position_static` or
            :class:`spatialtis.plotting.base.graph_position_interactive`

    {pyecharts_tips}

    """
    if community_key is None:
        community_key = ANALYSIS["cell_community"].last_used_key
    if centroid_key is None:
        centroid_key = CONFIG.CENTROID_KEY
    if neighbors_key is None:
        neighbors_key = CONFIG.NEIGHBORS_KEY

    df = data.obs.query("&".join([f"({k}=='{v}')" for k, v in roi.items()]))

    nodes_types = df[community_key].tolist()
    commus = []
    for commu, count in Counter(nodes_types).items():
        if count >= min_cells:
            commus.append(commu)

    df = df.reset_index(drop=True)
    xdf = df[df[community_key].isin(commus)]
    xdf = xdf.reset_index()

    need_eval_nodes = isinstance(xdf[centroid_key][0], str)
    need_eval_neighs = isinstance(xdf[neighbors_key][0], str)
    if need_eval_nodes:
        nodes = [literal_eval(n) for n in xdf[centroid_key]]
    else:
        nodes = [n for n in xdf[centroid_key]]
    if need_eval_neighs:
        neighs = [literal_eval(n) for n in xdf[neighbors_key]]
    else:
        neighs = [n for n in xdf[neighbors_key]]
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

    plot_options["saved_name"] = "community_map_" + ",".join(
        [f"{k}={v}" for k, v in roi.items()]
    )
    if use == "interactive":
        return graph_position_interactive(
            nodes, edges, edges_types=edges_types, **plot_options
        )
    else:
        return graph_position_static(
            nodes, edges, edges_types=edges_types, **plot_options
        )
