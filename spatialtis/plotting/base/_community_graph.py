from pathlib import Path
from typing import Optional, Sequence, Union

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import LineCollection
from matplotlib.lines import Line2D
from pyecharts import options as opts
from pyecharts.charts import Graph

from spatialtis.config import CONFIG
from spatialtis.utils import reuse_docstring

from .palette import get_colors
from .save import save_pyecharts


@reuse_docstring()
def graph_plot_interactive(
    nodes: Sequence,
    edges: Sequence,
    nodes_types: Optional[Sequence] = None,
    edges_types: Optional[Sequence] = None,
    node_size: Union[float, int] = 1,
    edge_size: Union[float, int] = 0.5,
    size: Sequence = (800, 800),
    renderer: str = "canvas",
    theme: str = "white",
    display: Optional[bool] = None,
    return_plot: bool = False,
    title: Optional[str] = None,
    save: Union[str, Path, None] = None,
):
    """(pyecharts) Graph visualization

    Args:
        nodes: List of nodes
        edges: List of edges
        nodes_types: The type of every node
        edges_types: The type of every edge
        node_size: The size of node
        edge_size: The size of edge width
        size: {size}
        renderer: {renderer}
        theme: {theme}
        title: {title}
        display: {display}
        save: {size}
        return_plot: {return_plot}

    """
    nodes_data = []
    edges_data = []
    categories = []

    if nodes_types is not None:
        n_unitypes = np.unique(nodes_types)
        for c in n_unitypes:
            categories.append(opts.GraphCategory(name=c, symbol="circle"))

    if edges_types is not None:
        e_unitypes = np.unique(edges_types)
        edges_colors = get_colors(len(e_unitypes), ["Set3"])
        edges_colormap = dict(zip(e_unitypes, edges_colors))

    for i, (x, y) in enumerate(nodes):

        node_config = dict(
            name=str(i),
            x=x,
            y=-y,
            label_opts=opts.LabelOpts(is_show=False),
            symbol_size=node_size,
        )

        if nodes_types is not None:
            node_config["category"] = nodes_types[i]

        nodes_data.append(opts.GraphNode(**node_config))

    for i, (source, target) in enumerate(edges):

        edge_config = dict(source=str(source), target=str(target),)

        if edges_types is not None:
            edge_config["linestyle_opts"] = opts.LineStyleOpts(
                width=edge_size, color=edges_colormap[edges_types[i]],
            )
        else:
            edge_config["linestyle_opts"] = opts.LineStyleOpts(width=edge_size,)

        edges_data.append(opts.GraphLink(**edge_config))

    g = Graph(
        init_opts=opts.InitOpts(
            width=f"{size[0]}px", height=f"{size[1]}px", renderer=renderer, theme=theme,
        )
    )
    g.add(
        "",
        nodes_data,
        edges_data,
        categories,
        layout="none",
        is_rotate_label=True,
        edge_label=opts.LabelOpts(is_show=False),
        tooltip_opts=opts.TooltipOpts(formatter="{c}"),
    ).set_global_opts(
        title_opts=opts.TitleOpts(title=title),
        # visualmap_opts=opts.VisualMapOpts(range_color=palette),
        legend_opts=opts.LegendOpts(
            type_="scroll",
            orient="vertical",
            pos_left="left",
            pos_top="center",
            item_width=5,
            legend_icon="circle",
            padding=3,
            textstyle_opts={"fontSize": 8},
        ),
        toolbox_opts=opts.ToolboxOpts(
            feature={
                "saveAsImage": {"title": "save", "pixelRatio": 5,},
                "brush": {},
                "restore": {},
            },
        ),
    )

    if save is not None:
        save_pyecharts(g, save)

    if display is None:
        if CONFIG.WORKING_ENV is None:
            display = False
        else:
            display = True
    if display:
        g.load_javascript()
        return g.render_notebook()

    if return_plot:
        return g


@reuse_docstring()
def graph_plot(
    nodes,
    edges,
    nodes_types=None,
    edges_types=None,
    node_size: Union[float, int] = 1,
    edge_size: Union[float, int] = 0.7,
    size: Sequence = (15, 15),
    palette: Optional[Sequence] = None,
    display: Optional[bool] = None,
    return_plot: bool = False,
    title: Optional[str] = None,
    save: Union[str, Path, None] = None,
):
    """(matplotlib) Graph visualization

    Args:
        nodes: List of nodes
        edges: List of edges
        nodes_types: The type of every node
        edges_types: The type of every edge
        node_size: The size of node
        edge_size: The size of edge width
        size: {size}
        palette: {palette}
        title: {title}
        display: {display}
        save: {save}
        return_plot: {return_plot}

    """
    if nodes_types is not None:
        n_unitypes = np.unique(nodes_types)
        nodes_colors = get_colors(len(n_unitypes), ["Category20"])
        nodes_colormap = dict(zip(n_unitypes, nodes_colors))

        nodes_legend = []
        for label in n_unitypes:
            nodes_legend.append(
                Line2D(
                    (),
                    (),
                    color="white",
                    marker="o",
                    markerfacecolor=nodes_colormap[label],
                    label=label,
                    markersize=8,
                ),
            )
    if edges_types is not None:
        e_unitypes = np.unique(edges_types)
        if palette is None:
            edges_colors = get_colors(len(e_unitypes), ["Spectral", "Set3"])
        else:
            edges_colors = get_colors(len(e_unitypes), palette)
        edges_colormap = dict(zip(e_unitypes, edges_colors))

    fig, ax = plt.subplots(figsize=size)

    line_cols = []
    line_colors = []
    for i, (e1, e2) in enumerate(edges):
        line_cols.append([nodes[e1], nodes[e2]])

        if edges_types is not None:
            line_color = edges_colormap[edges_types[i]]
        else:
            line_color = "#cccccc"
        line_colors.append(line_color)
    segs = LineCollection(line_cols, zorder=-1, linewidth=edge_size, colors=line_colors)

    ax.add_collection(segs)

    for i, (n1, n2) in enumerate(nodes):
        plot_config = dict(s=node_size)
        if nodes_types is not None:
            plot_config["c"] = nodes_colormap[nodes_types[i]]
        else:
            plot_config["c"] = "#CB1B45"
        ax.scatter(n1, n2, **plot_config)

    if nodes_types is not None:
        nlegend = ax.legend(
            handles=nodes_legend,
            loc="center left",
            bbox_to_anchor=(1, 0.5),
            frameon=False,
            labelspacing=1.2,
        )
        # nlegend._legend_box.align = "left"
        ax.add_artist(nlegend)

    plt.axis("off")

    if display is None:
        if CONFIG.WORKING_ENV is None:
            display = False
        else:
            display = True
    if not display:
        plt.close()

    if title:
        plt.title(title)

    if save:
        plt.savefig(save, dpi=300, bbox_inches="tight")

    if return_plot:
        return fig, ax
