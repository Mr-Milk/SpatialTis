from pathlib import Path
from typing import Mapping, Optional, Sequence, Union

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import PatchCollection
from matplotlib.lines import Line2D
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

from .palette import get_linear_colors


def dot_matrix(
    matrix: np.array,
    dot_color: np.array,
    dot_size: np.array,
    xlabels: Sequence,
    ylabels: Sequence,
    xlabel_rotation: int = 90,
    ylabel_rotation: int = 0,
    cbar_mapper: Mapping = None,
    cbar_legend_title: Optional[str] = None,
    size_legend_title: Optional[str] = None,
    color_legend_title: Optional[str] = None,
    color_legend_text: Optional[Sequence[str]] = None,
    palette: Optional[Sequence] = None,
    display: bool = True,
    return_plot: bool = False,
    title: Optional[str] = None,
    save: Union[str, Path, None] = None,
):
    """(matplotlib) Dot matrix plot

    Args:
        matrix: the value to control the color of matrix
        dot_color: the value to control the dot color
        dot_size: the value to control the dot size
        xlabels: the labels marked on x axis
        ylabels: the labels marked on y axis
        xlabel_rotation: rotate the x label
        ylabel_rotation: rotate the y label
        cbar_mapper: how to map the colorbar to annotation text
        cbar_legend_title: title of colorbar
        size_legend_title: title of size legend
        color_legend_title: title of color legend
        color_legend_text: the text in the color legend
        palette: the color
        display: whether to show the plot
        return_plot: whether to return the plot
        title: title of the plot
        save: the path to save your plot

    """
    if palette is not None:
        cmap = get_linear_colors(palette)
    else:
        cmap = "RdYlGn"

    N, M = np.asarray(matrix).shape
    x, y = np.meshgrid(np.arange(M), np.arange(N))

    s = np.asarray(dot_size)
    c = np.asarray(dot_color)

    fig, ax = plt.subplots(figsize=(M / 2, N / 2))

    R = s / s.max() / 2 * 0.9

    circles = [plt.Circle((j, i), radius=r,) for r, j, i in zip(R.flat, x.flat, y.flat)]
    circ_col = PatchCollection(circles, array=c.flatten(), cmap=cmap, alpha=0.5)

    rects = [plt.Rectangle((j - 0.5, i - 0.5), 1, 1,) for j, i in zip(x.flat, y.flat)]
    rect_col = PatchCollection(rects, array=matrix.flatten(), cmap=cmap, alpha=0.5)

    ax.add_collection(rect_col)
    ax.add_collection(circ_col)

    ax.set(
        xticks=np.arange(M),
        yticks=np.arange(N),
        xticklabels=xlabels,
        yticklabels=ylabels,
    )
    ax.set_xticks(np.arange(M + 1) - 0.5, minor=True)
    ax.set_yticks(np.arange(N + 1) - 0.5, minor=True)

    sm = dot_size.max()
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

    legend = ax.legend(
        handles=legends,
        title=size_legend_title,
        loc="upper left",
        bbox_to_anchor=(1.05, 0, 0, 1),
        frameon=False,
        markerscale=21,
        labelspacing=1.2,
    )
    legend._legend_box.align = "left"
    ax.add_artist(legend)

    colors = mpl.cm.get_cmap("RdYlGn")
    clegends = [
        Line2D(
            (),
            (),
            color="white",
            marker="o",
            label=color_legend_text[0],
            markerfacecolor=colors(1.0),
            markersize=0.7,
            alpha=0.5,
        ),
        Line2D(
            (),
            (),
            color="white",
            marker="o",
            label=color_legend_text[1],
            markerfacecolor=colors(0.5),
            markersize=0.7,
            alpha=0.5,
        ),
        Line2D(
            (),
            (),
            color="white",
            marker="o",
            label=color_legend_text[2],
            markerfacecolor=colors(0.0),
            markersize=0.7,
            alpha=0.5,
        ),
    ]

    clegend = ax.legend(
        handles=clegends,
        title=color_legend_title,
        loc="upper left",
        bbox_to_anchor=(1.05, -0.5, 0.5, 1),
        frameon=False,
        markerscale=20,
        labelspacing=1.2,
    )
    clegend._legend_box.align = "left"
    ax.add_artist(clegend)

    # set colorbar
    rect_col.set_clim([-1, 1])

    # turn off outer frame
    for spine in plt.gca().spines.values():
        spine.set_visible(False)

    # turn off minor ticks
    plt.minorticks_off()

    plt.xticks(rotation=xlabel_rotation)
    plt.yticks(rotation=ylabel_rotation)

    axins = inset_axes(
        ax,
        width="5%",  # width = 5% of parent_bbox width
        height="10%",  # height : 50%
        loc="upper left",
        bbox_to_anchor=(1.07, -0.9, 1, 1),
        bbox_transform=ax.transAxes,
    )

    cbar = fig.colorbar(rect_col, cax=axins, ticks=list(cbar_mapper.keys()))
    cbar.ax.set_title(cbar_legend_title)
    cbar.ax.set_yticklabels(list(cbar_mapper.values()))
    cbar.ax.yaxis.set_tick_params(length=0)  # hide ticks
    cbar.set_alpha(0.5)  # be consistent with color in plot

    if title:
        plt.title(title)

    if save:
        fig.savefig(save, dpi=300, bbox_inches="tight")

    if not display:
        plt.close()

    if return_plot:
        return fig, ax
