from typing import Optional, Sequence, Union
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import PatchCollection
from matplotlib.lines import Line2D
from mpl_toolkits.axes_grid1.inset_locator import inset_axes


def dot_matrix(
        matrix,
        dot_color,
        dot_size,
        xlabels,
        ylabels,
        cbar_mapper=None,
        cbar_legend_title=None,
        size_legend_title=None,
        color_legend_title=None,
        color_legend_text=None,
        palette: Optional[Sequence] = None,
        display: bool = True,
        return_plot: bool = False,
        title: Optional[str] = None,
        save: Union[str, Path, None] = None,
):
    N, M = np.asarray(matrix).shape
    x, y = np.meshgrid(np.arange(M), np.arange(N))

    s = np.asarray(dot_size)
    c = np.asarray(dot_color)

    fig, ax = plt.subplots(figsize=(M / 2, N / 2))

    R = s / s.max() / 2 * 0.9

    circles = [plt.Circle((j, i), radius=r, ) for r, j, i in zip(R.flat, x.flat, y.flat)]
    circ_col = PatchCollection(circles, array=c.flatten(), cmap="RdYlGn", alpha=0.5)

    rects = [plt.Rectangle((j - 0.5, i - 0.5), 1, 1, ) for j, i in zip(x.flat, y.flat)]
    rect_col = PatchCollection(rects, array=matrix.flatten(), cmap="RdYlGn", alpha=0.5)

    ax.add_collection(rect_col)
    ax.add_collection(circ_col)

    ax.set(xticks=np.arange(M), yticks=np.arange(N),
           xticklabels=xlabels, yticklabels=ylabels)
    ax.set_xticks(np.arange(M + 1) - 0.5, minor=True)
    ax.set_yticks(np.arange(N + 1) - 0.5, minor=True)

    sm = dot_size.max()
    legends = [
        Line2D((), (), color="white", marker='o', markerfacecolor="black", label=f"{int(sm)}", markersize=1),
        Line2D((), (), color="white", marker='o', markerfacecolor="black", label=f"{int(sm/4*3)}", markersize=0.75),
        Line2D((), (), color="white", marker='o', markerfacecolor="black", label=f"{int(sm/2)}", markersize=0.5),
        Line2D((), (), color="white", marker='o', markerfacecolor="none", markeredgecolor="black", label="0",
               markersize=0.25),
    ]

    legend = ax.legend(handles=legends, title=size_legend_title, loc='upper left', bbox_to_anchor=(1.05, 0, 0, 1),
                       frameon=False, markerscale=21, labelspacing=1.2)
    legend._legend_box.align = "left"
    ax.add_artist(legend)

    colors = mpl.cm.get_cmap('RdYlGn')
    clegends = [
        Line2D((), (), color="white", marker='o', label=color_legend_text[0], markerfacecolor=colors(1.0), markersize=0.7, alpha=0.5),
        Line2D((), (), color="white", marker='o', label=color_legend_text[1], markerfacecolor=colors(0.5), markersize=0.7, alpha=0.5),
        Line2D((), (), color="white", marker='o', label=color_legend_text[2], markerfacecolor=colors(0.0), markersize=0.7, alpha=0.5),
        ]

    clegend = ax.legend(handles=clegends, title=color_legend_title, loc='upper left', bbox_to_anchor=(1.05, -0.5, 0.5, 1),
                        frameon=False, markerscale=20, labelspacing=1.2)
    clegend._legend_box.align = "left"
    ax.add_artist(clegend)

    # set colorbar
    rect_col.set_clim([-1, 1])

    # turn off outer frame
    for spine in plt.gca().spines.values():
        spine.set_visible(False)

    # turn off minor ticks
    plt.minorticks_off()

    axins = inset_axes(ax,
                       width="5%",  # width = 5% of parent_bbox width
                       height="10%",  # height : 50%
                       loc='upper left',
                       bbox_to_anchor=(1.07, -0.9, 1, 1),
                       bbox_transform=ax.transAxes,
                       )

    cbar = fig.colorbar(rect_col, cax=axins, ticks=list(cbar_mapper.keys()))
    cbar.ax.set_title(cbar_legend_title)
    cbar.ax.set_yticklabels(list(cbar_mapper.values()))
    cbar.ax.yaxis.set_tick_params(length=0)  # hide ticks
    cbar.set_alpha(0.5)  # be consistent with color in plot

    plt.show()
