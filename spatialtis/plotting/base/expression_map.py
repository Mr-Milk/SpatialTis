from typing import List, Sequence, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pyecharts.options as opts
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from pyecharts.charts import Bar3D

from spatialtis.plotting.abc import MatplotlibMixin, PyechartsMixin
from spatialtis.plotting.base.palette import get_linear_colors
from spatialtis.typing import Number
from spatialtis.utils import doc


@doc
class expression_map_3d(PyechartsMixin):
    """Visualize marker expression in 3D, Pyecharts

    Args:
        points: The locations of each cells
        expressions: The expression level for each cell
        axis_size: The length of x,y,z axis
        **plot_options: {plot_options}

    """

    def __init__(
        self,
        points: List[Tuple[Number, Number]],
        expressions: List[Number],
        axis_size: tuple = (100, 100, 80),
        **plot_options,
    ):
        super().__init__(**plot_options)

        if self.palette is None:
            self.palette = ["RdYlBu"]
        default_color = get_linear_colors(self.palette)

        zdata = []
        for cell, exp in zip(points, expressions):
            zdata.append([cell[0], cell[1], exp])

        zrange = sorted(zdata, key=lambda k: k[2])

        self.plot = Bar3D(
            init_opts=opts.InitOpts(
                width=f"{self.size[0]}px",
                height=f"{self.size[1]}px",
                theme=self.theme,
                renderer=self.renderer,
            )
        )

        self.plot.add(
            series_name="",
            shading="color",
            data=zdata,
            xaxis3d_opts=opts.Axis3DOpts(type_="value"),
            yaxis3d_opts=opts.Axis3DOpts(type_="value"),
            zaxis3d_opts=opts.Axis3DOpts(type_="value"),
            grid3d_opts=opts.Grid3DOpts(
                width=axis_size[1], height=axis_size[2], depth=axis_size[0]
            ),
        ).set_global_opts(
            visualmap_opts=opts.VisualMapOpts(
                dimension=2,
                max_=zrange[-1][2],
                min_=zrange[0][2],
                range_color=default_color,
            ),
            tooltip_opts=opts.TooltipOpts(is_show=False),
            toolbox_opts=opts.ToolboxOpts(
                feature={
                    "saveAsImage": {"title": "Save", "pixelRatio": 5,},
                    "restore": {"title": "Restore"},
                },
            ),
        )
        self.set_up()


@doc
class expression_map_static(MatplotlibMixin):
    """Visualize marker expression, Matplotlib

    Args:
        points: The locations of each cells
        expressions: The expression level for each cell
        cell_size: The size of cell
        **plot_options: {plot_options}

    """

    def __init__(
        self,
        points: List[Tuple[Number, Number]],
        expressions: List[Number],
        cell_size: int = 5,
        **plot_options,
    ):
        super().__init__(**plot_options)

        if self.palette is None:
            self.palette = "magma"
        else:
            if isinstance(self.palette, Sequence):
                self.palette = self.palette[0]

        x, y = [], []
        for p in points:
            x.append(p[0])
            y.append(p[1])
        self.fig, self.ax = plt.subplots(figsize=self.size)
        ss = self.ax.scatter(x=x, y=y, c=expressions, s=cell_size, cmap=self.palette)
        # Remove the legend and add a colorbar
        self.ax.set_aspect("equal")
        plt.axis("off")
        axins = inset_axes(
            self.ax,
            width="3%",  # width = 5% of parent_bbox width
            height="50%",  # height : 50%
            loc="lower left",
            bbox_to_anchor=(1, 0.05, 1, 1),
            bbox_transform=self.ax.transAxes,
            borderpad=0,
        )
        ticks = np.linspace(min(expressions), max(expressions), num=6)
        cbar = self.fig.colorbar(ss, cax=axins, ticks=ticks)
        cbar.set_ticklabels([round(i, 1) for i in ticks])
        self.set_up()
