from pathlib import Path
from typing import Mapping, Optional, Sequence, Union

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import PatchCollection
from matplotlib.lines import Line2D
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

from spatialtis.config import CONFIG
from spatialtis.utils import reuse_docstring

from .palette import get_linear_colors


@reuse_docstring()
class DotMatrix:
    """(matplotlib) Dot matrix plot

    Args:
        matrix: The value to control the color of matrix
        dot_color: The value to control the dot color
        dot_size: The value to control the dot size
        xlabels: {xlabels}
        ylabels: {ylabels}
        xlabel_rotation: {xlabel_rotation}
        ylabel_rotation: {ylabel_rotation}
        show_ticks: {show_ticks}
        xaxis_title: {xaxis_title}
        yaxis_title: {yaxis_title}
        matrix_cbar_text: How to map the rect colorbar to annotation text
        dot_cbar_text: How to map the dot colorbar to annotation text
        matrix_cbar_title: The title of the rect colorbar
        dot_cbar_title: The title of dot colorbar
        size_legend_title: Title of the size legend
        matrix_palette: Control the color of rect
        dot_palette: Control the color of dot
        display: {display}
        title: {title}
        save: {save}

    """

    def __repr__(self):
        return ""

    def __init__(
        self,
        matrix: Optional[Sequence] = None,
        dot_color: Optional[Sequence] = None,
        dot_size: Optional[Sequence] = None,
        xlabels: Optional[Sequence] = None,
        ylabels: Optional[Sequence] = None,
        xlabel_rotation: int = 90,
        ylabel_rotation: int = 0,
        show_ticks: bool = True,
        xaxis_title: Optional[str] = None,
        yaxis_title: Optional[str] = None,
        matrix_cbar_text: Optional[Mapping] = None,
        dot_cbar_text: Optional[Mapping] = None,
        matrix_cbar_title: Optional[str] = None,
        dot_cbar_title: Optional[str] = None,
        size_legend_title: Optional[str] = None,
        matrix_palette: Optional[Sequence] = None,
        dot_palette: Optional[Sequence] = None,
        display: Optional[bool] = None,
        title: Optional[str] = None,
        save: Union[str, Path, None] = None,
    ):

        # global arguments
        self.M = 0
        self.N = 0

        # set palette
        if matrix_palette is not None:
            self.matrix_cmap = get_linear_colors(matrix_palette)
        else:
            self.matrix_cmap = "RdYlGn"

        if dot_palette is not None:
            self.dot_cmap = get_linear_colors(dot_palette)
        else:
            self.dot_cmap = "PiYG"

        if xlabels is not None:
            self.xlabels = xlabels
        else:
            raise ValueError("Please specific `xlabels`")

        if ylabels is not None:
            self.ylabels = ylabels
        else:
            self.ylabels = xlabels

        # save args
        if (matrix is None) & (dot_size is None):
            raise ValueError("Need either `matrix` or `dot_size`")

        if dot_size is not None:
            self.dot_size = np.asarray(dot_size)
            self.size_max = self.dot_size.flatten().max()
        if dot_color is not None:
            self.dot_color = np.asarray(dot_color)

        if matrix is not None:
            self.matrix = np.asarray(matrix)
            self.N, self.M = self.matrix.shape
        else:
            self.N, self.M = self.dot_size.shape

        self.matrix_cbar_text = matrix_cbar_text
        self.dot_cbar_text = dot_cbar_text
        self.matrix_cbar_title = matrix_cbar_title
        self.dot_cbar_title = dot_cbar_title
        self.size_legend_title = size_legend_title

        self.x, self.y = np.meshgrid(np.arange(self.M), np.arange(self.N))
        self.fig, self.ax = plt.subplots(figsize=(self.M / 2, self.N / 2))
        self.ax.set_aspect("equal")

        # turn off minor ticks
        # self.ax.minorticks_off()
        self.ax.set(
            xticks=np.arange(self.M),
            yticks=np.arange(self.N),
            xticklabels=self.xlabels,
            yticklabels=self.ylabels,
        )
        self.ax.set_xticks(np.arange(self.M + 1) - 0.5, minor=True)
        self.ax.set_yticks(np.arange(self.N + 1) - 0.5, minor=True)
        if not show_ticks:
            self.ax.xaxis.set_ticks_position("none")
            self.ax.yaxis.set_ticks_position("none")
        self.ax.set_xlabel(xaxis_title)
        self.ax.set_ylabel(yaxis_title)

        plt.xticks(rotation=xlabel_rotation)
        plt.yticks(rotation=ylabel_rotation)

        if title:
            plt.title(title)

        # turn off outer frame
        for spine in self.ax.spines.values():
            spine.set_visible(False)

        if matrix is not None:
            self._draw_rect()
        if dot_size is not None:
            self._draw_dot()

        if save:
            self.fig.savefig(save, dpi=300, bbox_inches="tight")

        if display is None:
            if CONFIG.WORKING_ENV is None:
                display = False
            else:
                display = True
        if not display:
            plt.close()

    def _draw_rect(self):

        rects = [
            plt.Rectangle((j - 0.5, i - 0.5), 1, 1,)
            for j, i in zip(self.x.flat, self.y.flat)
        ]
        color_array = self.matrix.flatten()
        cmax = np.max(color_array)
        cmin = np.min(color_array)
        rect_col = PatchCollection(
            rects, array=color_array, cmap=self.matrix_cmap, alpha=0.5
        )
        rect_col.set_clim(cmin, cmax)
        self.ax.add_collection(rect_col)
        self._draw_cbar(
            rect_col,
            list(np.linspace(cmin, cmax, num=len(self.dot_cbar_text))),
            self.matrix_cbar_text,
            (1.07, -0.9, 1, 1),
            title=self.matrix_cbar_title,
        )

    def _draw_dot(self):

        R = self.dot_size / self.size_max / 2 * 0.9
        circles = [
            plt.Circle((j, i), radius=r,)
            for r, j, i in zip(R.flat, self.x.flat, self.y.flat)
        ]
        color_array = self.dot_color.flatten()
        cmax = np.max(color_array)
        cmin = np.min(color_array)
        circ_col = PatchCollection(
            circles, array=color_array, cmap=self.dot_cmap, alpha=0.5
        )
        circ_col.set_clim(cmin, cmax)
        self.ax.add_collection(circ_col)
        self._draw_dot_size_legend()
        self._draw_cbar(
            circ_col,
            list(np.linspace(cmin, cmax, num=len(self.dot_cbar_text))),
            self.dot_cbar_text,
            (1.07, -0.5, 1, 1),
            title=self.dot_cbar_title,
        )

    def _draw_dot_size_legend(self):
        sm = self.size_max

        legends = [
            Line2D(
                (),
                (),
                color="white",
                marker="o",
                markerfacecolor="black",
                label=f"{int(sm)}",
                markersize=1,
            ),
            Line2D(
                (),
                (),
                color="white",
                marker="o",
                markerfacecolor="black",
                label=f"{int(sm / 4 * 3)}",
                markersize=0.75,
            ),
            Line2D(
                (),
                (),
                color="white",
                marker="o",
                markerfacecolor="black",
                label=f"{int(sm / 2)}",
                markersize=0.5,
            ),
            Line2D(
                (),
                (),
                color="white",
                marker="o",
                markerfacecolor="none",
                markeredgecolor="black",
                label="0",
                markersize=0.25,
            ),
        ]

        legend = self.ax.legend(
            handles=legends,
            title=self.size_legend_title,
            loc="upper left",
            bbox_to_anchor=(1.05, 0, 0, 1),
            frameon=False,
            markerscale=21,
            labelspacing=1.2,
        )

        legend._legend_box.align = "left"
        self.ax.add_artist(legend)

    def _draw_cbar(self, collections, ticks, ticks_text, position, title=None):
        axins = inset_axes(
            self.ax,
            width="5%",  # width = 5% of parent_bbox width
            height="10%",  # height : 50%
            loc="upper left",
            bbox_to_anchor=position,
            bbox_transform=self.ax.transAxes,
        )
        cbar = self.fig.colorbar(collections, cax=axins, ticks=ticks)
        cbar.ax.set_title(title)
        cbar.ax.set_yticklabels(list(ticks_text))
        cbar.ax.yaxis.set_tick_params(length=0)  # hide ticks
        cbar.set_alpha(0.5)  # be consistent with color in plot
