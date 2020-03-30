from typing import Optional, Union, Sequence, Mapping
from collections import Counter
from pathlib import Path

import pandas as pd
import numpy as np

from pyecharts import options as opts
from pyecharts.charts import Graph

from .palette import get_linear_colors
from ._save import save_pyecharts


# mapper = {-1: "Avoidance", 1: "Association"}
def cc_interactions(df: pd.DataFrame,
                    mapper: Mapping,
                    order: Optional[Sequence] = None,
                    repulsion: Union[float, int] = 80000,
                    gravity: Union[float, int] = 0.2,
                    size: Sequence = (800, 800),
                    renderer: str = 'canvas',
                    theme: str = 'white',
                    edges_colors: Optional[Sequence] = None,
                    display: bool = True,
                    return_plot: bool = False,
                    title: Optional[str] = None,
                    save: Union[str, Path, None] = None,
                    ):
    if edges_colors is None:
        edges_colors = ["RdBu"]

    nodes_data = []
    edges_data = []

    edges_colors_range = get_linear_colors(edges_colors)
    edges_width_range = np.arange(1, 101) / 10
    nodes_size_range = np.arange(1, 5) * 10

    nodes = np.unique(list(df.columns.to_numpy()))

    edges = np.unique(df.to_numpy())
    if order is None:
        order = edges

    graphs_nodes = {n: 0 for n in nodes}
    graphs_edges = {}
    # (source, target): (type, value, sum) | (interaction_type, color, line_width)

    for i, (label, data) in enumerate(df.items()):
        c_data = sorted([(k, v) for k, v in Counter(data).items() if k in order], key=lambda x: x[1], reverse=True)
        graphs_edges[label] = (c_data[0][0], c_data[0][1], sum([i[1] for i in c_data]))

    # the bigger the nodes, the cell tend to interact with more cells
    for e, v in graphs_edges.items():
        if v[0] == 1:
            graphs_nodes[e[1]] += v[-1]

    # nodes values
    nv = sorted(list(graphs_nodes.values()))
    nodes_size_normf = (nv[-1] - nv[0]) / (nodes_size_range[-1] - nodes_size_range[0])
    nodes_size = {k: ((v - nv[0]) / nodes_size_normf + nodes_size_range[0]) for k, v in graphs_nodes.items()}

    # link width
    lw = sorted([v[-1] for v in graphs_edges.values()])
    lw_normf = (lw[-1] - lw[0]) / (edges_width_range[-1] - edges_width_range[0])
    link_width = {k: ((v[-1] - lw[0]) / lw_normf + edges_width_range[0]) for k, v in graphs_edges.items()}

    # link color
    lc = sorted([v[1] for v in graphs_edges.values()])

    lc_normf = (lc[-1] - lc[0]) / (len(edges_colors_range))

    link_color = {k: edges_colors_range[int((v[1] - lc[0]) // lc_normf)] for k, v in graphs_edges.items()}

    nodes_data = []
    edges_data = []

    for n, v in graphs_nodes.items():
        nodes_data.append(
            opts.GraphNode(name=n,
                           symbol_size=nodes_size[n],
                           value=f"{n}:{v}",
                           )
        )

    for e, v in graphs_edges.items():
        edges_data.append(
            opts.GraphLink(source=e[1],
                           target=e[0],
                           symbol=['', 'arrow'],
                           value=f"{mapper[v[0]]} {v[1]}/{v[-1]}",
                           symbol_size=link_width[e] * 2,
                           linestyle_opts=opts.LineStyleOpts(width=link_width[e], curve=0.2, color=link_color[e]),
                           label_opts=opts.LabelOpts(is_show=False, )
                           )
        )

    g = Graph(init_opts=opts.InitOpts(width=f"{size[0]}px",
                                      height=f"{size[1]}px",
                                      renderer=renderer,
                                      theme=theme,
                                      ))
    g.add("",
          nodes_data,
          edges_data,
          layout="force",
          repulsion=repulsion,
          gravity=gravity,
          is_rotate_label=True,
          edge_label=opts.LabelOpts(is_show=False),
          tooltip_opts=opts.TooltipOpts(formatter="{c}"),
          ).set_global_opts(
        title_opts=opts.TitleOpts(title=title),
        toolbox_opts=opts.ToolboxOpts(feature={"saveAsImage": {"title": "save", "pixelRatio": 5, }, }, )
    )

    if save is not None:
        save_pyecharts(g, save)

    if display:
        g.load_javascript()
        return g.render_notebook()

    if return_plot:
        return g
