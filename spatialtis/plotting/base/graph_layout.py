from typing import List, Optional, Tuple

import matplotlib as mpl
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from pyecharts import options as opts
from pyecharts.charts import Graph

from spatialtis.plotting.abc import MatplotlibMixin, PyechartsMixin
from spatialtis.plotting.base.palette import get_linear_colors
from spatialtis.typing import Number
from spatialtis.utils import doc


@doc
class graph_layout_interactive(PyechartsMixin):
    """Graph with layout

    Args:
        edges: [source, target, width, color]
        max_width: The max width of the edges
        min_width: The min width of the edges
        curve: The curveness of edges
        width_vmin: Minimum for edge width scaling
        width_vmax: Maximum for edge width scaling
        color_vmin: Minimum for edge color scaling
        color_vmax: Maximum for edge color scaling
        cbar_text: Text to put on colorbar
        layout: "force" or "circular" (Default: "force")
        directed: If true, the edges will have arrow,
        **plot_options: {plot_options}

    """

    def __init__(
        self,
        edges: List[Tuple[str, str, Number, Number]],
        max_width: int = 5,
        min_width: int = 0,
        curve: float = 0.2,
        width_vmin: Optional[int] = None,
        width_vmax: Optional[int] = None,
        color_vmin: Optional[int] = None,
        color_vmax: Optional[int] = None,
        cbar_text: Optional[List] = None,
        layout: str = "force",
        directed: bool = False,
        **plot_options,
    ):
        super().__init__(**plot_options)
        plot_nodes = []
        plot_edges = []

        nodes_data = []
        edges_data = []
        edges_width = []
        edges_colors = []

        for e in edges:
            if e[0] != e[1]:
                nodes_data.append(e[0])
                nodes_data.append(e[1])
                edges_data.append((str(e[0]), str(e[1])))
                edges_width.append(e[2])
                edges_colors.append(e[3])

        nodes_data = pd.unique(nodes_data)
        for n in nodes_data:
            plot_nodes.append(opts.GraphNode(name=n))

        if self.palette is None:
            self.palette = ["RdBu"]
        color_pool = get_linear_colors(self.palette)

        if width_vmin is None:
            width_vmin = min(edges_width)
        if width_vmax is None:
            width_vmax = max(edges_width)
        if color_vmin is None:
            color_vmin = min(edges_colors)
        if color_vmax is None:
            color_vmax = max(edges_colors)

        if width_vmin == width_vmax:
            edges_width = [max_width for _ in range(len(edges_data))]
        else:
            edges_width.append(width_vmin)
            edges_width.append(width_vmax)
            edges_width = pd.cut(
                edges_width,
                10 * (max_width - min_width),
                labels=np.array(range(10 * (max_width - min_width))) / 10,
            ).tolist()
            edges_width.pop()
            edges_width.pop()

        if color_vmin == color_vmax:
            edges_colors = [color_pool[0] for _ in range(len(edges_data))]
        else:
            edges_colors.append(color_vmin)
            edges_colors.append(color_vmax)
            edges_colors = pd.cut(
                edges_colors, len(color_pool), labels=color_pool
            ).tolist()
            edges_colors.pop()
            edges_colors.pop()

        symbol = "arrow" if directed else None
        for e, w, c in zip(edges_data, edges_width, edges_colors):
            plot_edges.append(
                opts.GraphLink(
                    source=e[0],
                    target=e[1],
                    symbol=symbol,
                    linestyle_opts=opts.LineStyleOpts(width=w, color=c, curve=curve),
                )
            )

        self.plot = Graph(
            init_opts=opts.InitOpts(
                width=f"{self.size[0]}px",
                height=f"{self.size[1]}px",
                renderer=self.renderer,
                theme=self.theme,
            )
        )
        self.plot.add(
            "", plot_nodes, plot_edges, layout=layout, gravity=0.1,
        ).set_global_opts(
            visualmap_opts=opts.VisualMapOpts(
                range_color=color_pool, is_calculable=False, range_text=cbar_text,
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
class graph_layout_static(MatplotlibMixin):
    """Graph with layout, Matplotlib[NetworkX]

    Args:
        edges: [source, target, width, color]
        max_width: The max width of the edges
        min_width: The min width of the edges
        curve: The curveness of edges
        width_vmin: Minimum for edge width scaling
        width_vmax: Maximum for edge width scaling
        color_vmin: Minimum for edge color scaling
        color_vmax: Maximum for edge color scaling
        cbar_text: Text to put on colorbar
        layout: Methods from
            `networkx.layout <https://networkx.org/documentation/stable/reference/drawing.html#module-networkx.drawing.layout>`_,
        directed: If true, the edges will have arrow,
        **plot_options: {plot_options}

    """

    def __init__(
        self,
        edges: List[Tuple[str, str, Number, Number]],
        max_width: int = 5,
        min_width: int = 0,
        curve: float = 0.2,
        width_vmin: Optional[int] = None,
        width_vmax: Optional[int] = None,
        color_vmin: Optional[int] = None,
        color_vmax: Optional[int] = None,
        cbar_text: Optional[List] = None,
        layout: str = "spring_layout",
        directed: bool = False,
        **plot_options,
    ):
        super().__init__(**plot_options)
        edges_data = []
        edges_width = []
        edges_colors = []

        for e in edges:
            if e[0] != e[1]:
                edges_data.append((e[0], e[1]))
                edges_width.append(e[2])
                edges_colors.append(e[3])

        if width_vmin is None:
            width_vmin = min(edges_width)
        if width_vmax is None:
            width_vmax = max(edges_width)
        if width_vmin == width_vmax:
            edges_width = max_width
        else:
            edges_width.append(width_vmin)
            edges_width.append(width_vmax)
            edges_width = pd.cut(
                edges_width,
                10 * (max_width - min_width),
                labels=np.array(range(10 * (max_width - min_width))) / 10,
            ).tolist()
            edges_width.pop()
            edges_width.pop()

        edges_pool = {
            es: [w, c] for es, w, c in zip(edges_data, edges_width, edges_colors)
        }

        G = nx.DiGraph()
        G.add_edges_from(edges_data)
        if layout == "bipartite_layout":
            pos = getattr(nx.layout, layout).__call__(G, [e[0] for e in edges_data])
        else:
            pos = getattr(nx.layout, layout).__call__(G)

        self.fig, self.ax = plt.subplots(figsize=self.size)
        nx.draw_networkx_nodes(G, pos)
        symbol = "->" if directed else "-"
        if self.palette is None:
            self.palette = ["RdBu"]
        color_pool = get_linear_colors(self.palette)
        cmap = mpl.colors.LinearSegmentedColormap.from_list("", color_pool)
        edges_width = [edges_pool[e][0] for e in G.edges()]
        edges_colors = [edges_pool[e][1] for e in G.edges()]
        draw_edges = nx.draw_networkx_edges(
            G,
            pos,
            edge_color=edges_colors,
            edge_cmap=cmap,
            edge_vmin=color_vmin,
            edge_vmax=color_vmax,
            width=edges_width,
            arrowstyle=symbol,
            connectionstyle=f"arc3,rad={curve}",
        )
        nx.draw_networkx_labels(G, pos)
        axins = inset_axes(
            self.ax,
            width="5%",  # width = 5% of parent_bbox width
            height="10%",  # height : 50%
            loc="lower left",
            bbox_to_anchor=(1, 0.05, 1, 5),
            bbox_transform=self.ax.transAxes,
        )
        pc = mpl.collections.PatchCollection(draw_edges, cmap=cmap)
        edges_colors.append(color_vmin)
        edges_colors.append(color_vmax)
        pc.set_array(edges_colors)
        if cbar_text is not None:
            cbar = self.fig.colorbar(
                pc,
                cax=axins,
                ticks=np.linspace(color_vmin, color_vmax, num=len(cbar_text)),
            )
            cbar.set_ticklabels(cbar_text)
        else:
            ticks = np.linspace(color_vmin, color_vmax, num=4)
            cbar = self.fig.colorbar(pc, cax=axins, ticks=ticks)
            cbar.set_ticklabels([round(i, 1) for i in ticks])
        cbar.ax.set_title(self.legend_title)
        cbar.ax.yaxis.set_tick_params(length=0)  # hide ticks
        # cbar.set_alpha(0.5)  # be consistent with color in plot
        self.ax.grid(False)
        self.set_up()
