from typing import Optional

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import LineCollection
from matplotlib.lines import Line2D
from pyecharts import options as opts
from pyecharts.charts import Graph

from ...typing import Array, Number
from ...utils import doc
from ..abc import MatplotlibMixin, PyechartsMixin
from .palette import get_colors


@doc
class graph_position_interactive(PyechartsMixin):
    """Graph with fixed position

    Args:
        nodes: List of X, Y position of nodes
        edges: List of edges
        nodes_types: The type of every node
        edges_types: The type of every edge
        node_size: The size of node
        edge_size: The size of edge width
        **plot_options: {plot_options}

    """

    def __init__(
        self,
        nodes: Array,
        edges: Array,
        nodes_types: Optional[Array] = None,
        edges_types: Optional[Array] = None,
        node_size: Number = 1,
        edge_size: Number = 0.5,
        **plot_options,
    ):
        super().__init__(**plot_options)
        nodes_data = []
        edges_data = []
        categories = []

        if nodes_types is not None:
            n_unitypes = np.unique(nodes_types)
            for c in n_unitypes:
                categories.append(opts.GraphCategory(name=c, symbol="circle"))

        if self.palette is None:
            self.palette = ["Category20"]

        if edges_types is not None:
            e_unitypes = np.unique(edges_types)
            edges_colors = get_colors(len(e_unitypes), self.palette)
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

        self.plot = Graph(
            init_opts=opts.InitOpts(
                width=f"{self.size[0]}px",
                height=f"{self.size[1]}px",
                renderer=self.renderer,
                theme=self.theme,
            )
        )
        self.plot.add(
            "",
            nodes_data,
            edges_data,
            categories,
            layout="none",
            is_rotate_label=True,
            edge_label=opts.LabelOpts(is_show=False),
            tooltip_opts=opts.TooltipOpts(formatter="{c}"),
        ).set_global_opts(
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
                    "saveAsImage": {"title": "Save", "pixelRatio": 5,},
                    "restore": {"title": "Restore"},
                },
            ),
        )
        self.set_up()


@doc
class graph_position_static(MatplotlibMixin):
    """Graph with fixed position, Matplotlib[NetworkX]

    Args:
        nodes: List of nodes
        edges: List of edges
        nodes_types: The type of every node
        edges_types: The type of every edge
        node_size: The size of node
        edge_size: The size of edge width
        **plot_options: {plot_options}

    """

    def __init__(
        self,
        nodes: Array,
        edges: Array,
        nodes_types: Optional[Array] = None,
        edges_types: Optional[Array] = None,
        node_size: Number = 1,
        edge_size: Number = 0.5,
        **plot_options,
    ):
        super().__init__(**plot_options)
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
            if self.palette is None:
                edges_colors = get_colors(len(e_unitypes), ["Spectral", "Set3"])
            else:
                edges_colors = get_colors(len(e_unitypes), self.palette)
            edges_colormap = dict(zip(e_unitypes, edges_colors))

        self.fig, self.ax = plt.subplots(figsize=self.size)

        line_cols = []
        line_colors = []
        for i, (e1, e2) in enumerate(edges):
            line_cols.append([nodes[e1], nodes[e2]])

            if edges_types is not None:
                line_color = edges_colormap[edges_types[i]]
            else:
                line_color = "#cccccc"
            line_colors.append(line_color)
        segs = LineCollection(
            line_cols, zorder=-1, linewidth=edge_size, colors=line_colors
        )

        self.ax.add_collection(segs)

        for i, (n1, n2) in enumerate(nodes):
            plot_config = dict(s=node_size)
            if nodes_types is not None:
                plot_config["c"] = nodes_colormap[nodes_types[i]]
            else:
                plot_config["c"] = "#CB1B45"
            self.ax.scatter(n1, n2, **plot_config)

        if nodes_types is not None:
            nlegend = self.ax.legend(
                title=self.legend_title,
                handles=nodes_legend,
                ncol=2,
                loc="lower left",
                bbox_to_anchor=(1, 0.05),
                frameon=False,
                borderaxespad=0,
            )
            # nlegend._legend_box.align = "left"
            self.ax.add_artist(nlegend)
        plt.axis("off")
        self.set_up()
