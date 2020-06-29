from collections import Counter
from pathlib import Path
from typing import Mapping, Optional, Sequence, Union

import numpy as np
import pandas as pd
from pyecharts import options as opts
from pyecharts.charts import Graph

from ._save import save_pyecharts
from .palette import get_linear_colors


# mapper = {-1: "Avoidance", 1: "Association"}
# the direction of interactions is from first element to second element in tuple
def cc_interactions(
    df: pd.DataFrame,
    mapper: Mapping,
    repulsion: Union[float, int] = 80000,
    gravity: Union[float, int] = 0.2,
    threshold: float = 0.5,
    layout: str = "circular",
    renderer: str = "canvas",
    theme: str = "light",
    edges_colors: Optional[Sequence] = None,
    size: Sequence = (800, 800),
    title: Optional[str] = None,
    display: bool = True,
    save: Union[str, Path, None] = None,
    return_plot: bool = False,
):
    """(pyecharts) wrapper function for cell-cell interaction in graph

    Args:
        df: input data
        mapper: how to map your data to certain string
        repulsion: the force that push the node away
        gravity: the force that draw the node close
        threshold: the threhold to control vertice to display
        layout: "circular" or "force"
        renderer: "canvas" or "svg"
        theme: https://pyecharts.org/#/zh-cn/themes
        edges_colors:
        size: size of plot in pixels
        title: title of the plot
        display: whether to display the plot
        save: the path to save your plot
        return_plot: whether to return the plot instance

    Returns:

    """
    if edges_colors is None:
        edges_colors = ["RdBu"]

    edges_colors_range = get_linear_colors(edges_colors)
    edges_width_range = np.arange(1, 20) / 10
    nodes_size_range = np.arange(1, 5) * 10

    nodes = np.unique(list(df.columns.to_numpy()))

    edges = np.unique(df.to_numpy())

    graphs_nodes = {n: 0 for n in nodes}
    edges_data = {}
    graphs_edges = {}
    # (source, target): (associate or avoidance, % of ROIs that are sign in all ROIs.)

    for i, (label, data) in enumerate(df.items()):
        data = list(data.to_numpy())
        if label not in graphs_nodes.keys():
            reverse_label = (label[1], label[0])
            if reverse_label in graphs_nodes.keys():
                edges_data[reverse_label] += data
            else:
                edges_data[label] = data

    for c, data in edges_data.items():
        c_data = sorted(
            [(k, v) for k, v in Counter(data).items() if k in edges],
            key=lambda x: x[1],
            reverse=True,
        )
        proportions = c_data[0][1] / len(data)
        # if no interaction, drop the edge
        if (c_data[0][0] != 0) & (proportions > threshold):
            graphs_edges[c] = (c_data[0][0], proportions, c_data[0][1])

    # the bigger the nodes, the cell tend to interact with more cells in more ROIs
    for e, v in graphs_edges.items():
        if v[0] == 1:
            graphs_nodes[e[1]] += v[-1]

    # nodes values
    nv = sorted(list(graphs_nodes.values()))
    nodes_size_normf = (nv[-1] - nv[0]) / (nodes_size_range[-1] - nodes_size_range[0])
    nodes_size = {
        k: ((v - nv[0]) / nodes_size_normf + nodes_size_range[0])
        for k, v in graphs_nodes.items()
    }

    # link width
    """
    lw = sorted([v[-1] for v in graphs_edges.values()])
    lw_normf = (lw[-1] - lw[0]) / (edges_width_range[-1] - edges_width_range[0])
    link_width = {k: ((v[-1] - lw[0]) / lw_normf + edges_width_range[0]) for k, v in graphs_edges.items()}
    """

    # link color
    lc = len(edges_colors_range)
    link_color = {}

    for k, v in graphs_edges.items():
        if v[0] == 1:
            link_color[k] = edges_colors_range[int(lc * v[1]) - 1]
        elif v[0] == -1:
            link_color[k] = edges_colors_range[int(lc * (1 - v[1])) - 1]

    nodes_data = []
    edges_data = []

    for n, v in graphs_nodes.items():
        nodes_data.append(
            opts.GraphNode(name=n, symbol_size=nodes_size[n], value=f"{n}:{v}",)
        )

    for e, v in graphs_edges.items():
        edges_data.append(
            opts.GraphLink(
                source=e[0],
                target=e[1],
                # symbol=['', 'arrow'],
                value=f"{mapper[v[0]]}:{round(v[1], 2)} {e[0]}-{e[1]}",
                # symbol_size=link_width[e] * 2,
                linestyle_opts=opts.LineStyleOpts(
                    width=3, curve=0.2, color=link_color[e]
                ),
                label_opts=opts.LabelOpts(is_show=False,),
            )
        )

    g = Graph(
        init_opts=opts.InitOpts(
            width=f"{size[0]}px",
            height=f"{size[1]}px",
            renderer=renderer,
            theme=theme,
            animation_opts=opts.AnimationOpts(animation=False),
        )
    )
    g.add(
        "",
        nodes_data,
        edges_data,
        layout=layout,
        repulsion=repulsion,
        gravity=gravity,
        is_rotate_label=True,
        edge_label=opts.LabelOpts(is_show=False),
        tooltip_opts=opts.TooltipOpts(formatter="{c}"),
    ).set_global_opts(
        title_opts=opts.TitleOpts(title=title),
        toolbox_opts=opts.ToolboxOpts(
            feature={"saveAsImage": {"title": "save", "pixelRatio": 5,},},
        ),
    )

    if save is not None:
        save_pyecharts(g, save)

    if display:
        g.load_javascript()
        return g.render_notebook()

    if return_plot:
        return g
