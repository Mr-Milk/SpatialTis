from pathlib import Path
from typing import Mapping, Optional, Sequence, Union

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import PatchCollection
from matplotlib.lines import Line2D
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

from ..config import CONFIG
from .palette import get_linear_colors

# the relationship between label and arr should look like this

"""
    A B C D
A   1 2 3 4
B   1 2 3
C   1 2
D   1
"""


class TriDotMatrix:
    """(matplotlib) Triangular dot plot

    This input diagnol array should look like this

    .. code-block:: none

               A B C D
            A  1 2 3 4
            B  1 2 3
            C  1 2
            D  1

    Args:
        diagonal_block: the diagonal array for block
        diagonal_dot: the diagonal array for block
        labels: labels will be marked in x and y axis
        annotate: whether to show value on the dot
        xlabel_rotation: rotate x label
        ylabel_rotation: rotate y label
        legend_title: the title of legend
        color: color of the dot
        display: whether to show the plot
        title: title of the plot
        save: the path to save the plot

    """

    @staticmethod
    def make_flat(arr):
        flat = []
        for i in arr:
            for t in i:
                flat.append(t)
        return flat

    def __repr__(self):
        return ""

    def __init__(
        self,
        diagonal_block: Optional[Sequence] = None,
        diagonal_dot_color: Optional[Sequence] = None,
        diagonal_dot_size: Optional[Sequence] = None,
        labels: Optional[Sequence] = None,
        xlabel_rotation: int = 90,
        ylabel_rotation: int = 0,
        xaxis_title: Optional[str] = None,
        yaxis_title: Optional[str] = None,
        block_cbar_mapper: Optional[Mapping] = None,
        dot_cbar_mapper: Optional[Mapping] = None,
        block_cbar_title: Optional[str] = None,
        dot_cbar_title: Optional[str] = None,
        size_legend_title: Optional[str] = None,
        block_palette: Optional[Sequence] = None,
        dot_palette: Optional[Sequence] = None,
        display: Optional[bool] = None,
        title: Optional[str] = None,
        save: Union[str, Path, None] = None,
    ):

        # global arguments
        self.labels = labels
        self.M = self.N = len(labels)

        self.xlabel_rotation = xlabel_rotation
        self.ylabel_rotation = ylabel_rotation
        self.block_cbar_mapper = block_cbar_mapper
        self.dot_cbar_mapper = dot_cbar_mapper
        self.block_cbar_title = block_cbar_title
        self.dot_cbar_title = dot_cbar_title
        self.size_legend_title = size_legend_title

        # set palette
        if block_palette is not None:
            self.matrix_cmap = get_linear_colors(block_palette)
        else:
            self.matrix_cmap = "RdYlGn"

        if dot_palette is not None:
            self.dot_cmap = get_linear_colors(dot_palette)
        else:
            self.dot_cmap = "PiYG"

        self.block = diagonal_block
        self.dot_size = diagonal_dot_size
        self.dot_color = diagonal_dot_color

        flat_size = self.make_flat(diagonal_dot_size)
        self.circ_max = max(flat_size)

        self.fig, self.ax = plt.subplots(figsize=(self.N / 2, self.M / 2))
        self.ax.set_aspect("equal")
        self.ax.set(
            yticks=np.arange(self.M),
            xticks=np.arange(self.N),
            xticklabels=self.labels,
            yticklabels=self.labels[::-1],
        )
        self.ax.set_xticks(np.arange(self.M + 1) - 0.5, minor=True)
        self.ax.set_yticks(np.arange(self.N + 1) - 0.5, minor=True)
        self.ax.set_xlabel(xaxis_title)
        self.ax.set_ylabel(yaxis_title)

        plt.xticks(rotation=self.xlabel_rotation)
        plt.yticks(rotation=self.ylabel_rotation)

        if title:
            plt.title(title)

        # ax.xaxis.set_ticks_position("none")
        # ax.yaxis.set_ticks_position("none")

        for spine in plt.gca().spines.values():
            spine.set_visible(False)

        if diagonal_block is not None:
            self._draw_rect()
        if diagonal_dot_size is not None:
            self._draw_dot()

        if save:
            plt.savefig(save, dpi=300, bbox_inches="tight")

        if display is None:
            if CONFIG.WORKING_ENV is None:
                display = False
            else:
                display = True
        if not display:
            plt.close()

    def _draw_rect(self):

        flat_block = self.make_flat(self.block)
        rects = list()
        for i, arr in enumerate(self.block):
            for t, v in enumerate(arr[::-1]):
                rects.append(plt.Rectangle((i - 0.5, t - 0.5), 1, 1))
        rects_col = PatchCollection(
            rects, alpha=0.5, array=np.asarray(flat_block), cmap=self.matrix_cmap
        )
        self.ax.add_collection(rects_col)
        self._draw_cbar(
            rects_col,
            self.block_cbar_mapper,
            (1.07, -0.9, 1, 1),
            title=self.block_cbar_title,
        )

    def _draw_dot(self):

        flat_color = self.make_flat(self.dot_color)
        circles = list()
        for i, arr in enumerate(self.dot_size):
            for t, v in enumerate(arr[::-1]):
                if v != 0:
                    circles.append(
                        plt.Circle((i, t), radius=v / self.circ_max / 2 * 0.9)
                    )
        circ_col = PatchCollection(
            circles, alpha=0.5, array=np.asarray(flat_color), cmap=self.dot_cmap
        )
        self.ax.add_collection(circ_col)
        self._draw_dot_size_legend()
        if self.dot_cbar_mapper is not None:
            self._draw_cbar(
                circ_col,
                self.dot_cbar_mapper,
                (1.07, -0.5, 1, 1),
                title=self.dot_cbar_title,
            )

    def _draw_dot_size_legend(self):
        # size legends
        size_legends = [
            Line2D(
                (),
                (),
                color="white",
                marker="o",
                markerfacecolor="black",
                label=f"{int(self.circ_max)}",
                markersize=1.0,
            ),
            Line2D(
                (),
                (),
                color="white",
                marker="o",
                markerfacecolor="black",
                label=f"{int(self.circ_max / 4 * 3)}",
                markersize=0.75,
            ),
            Line2D(
                (),
                (),
                color="white",
                marker="o",
                markerfacecolor="black",
                label=f"{int(self.circ_max / 2)}",
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

        if self.size_legend_title is None:
            self.size_legend_title = "Size"

        self.ax.add_artist(
            self.ax.legend(
                handles=size_legends,
                loc="upper left",
                bbox_to_anchor=(1.07, -0.1, 0, 1),
                title=self.size_legend_title,
                markerscale=22,
                labelspacing=1.2,
                frameon=False,
            )
        )

    def _draw_cbar(self, collections, mapper, position, title=None):
        collections.set_clim([-1, 1])
        axins = inset_axes(
            self.ax,
            width="5%",  # width = 5% of parent_bbox width
            height="10%",  # height : 50%
            loc="upper left",
            bbox_to_anchor=position,
            bbox_transform=self.ax.transAxes,
        )
        cbar = self.fig.colorbar(collections, cax=axins, ticks=list(mapper.keys()))
        cbar.ax.set_title(title)
        cbar.ax.set_yticklabels(list(mapper.values()))
        cbar.ax.yaxis.set_tick_params(length=0)  # hide ticks
        cbar.set_alpha(0.5)  # be consistent with color in plot
