from pathlib import Path
from typing import Optional, Sequence, Union

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import PatchCollection
from matplotlib.lines import Line2D

from spatialtis.config import CONFIG
from spatialtis.utils import reuse_docstring


@reuse_docstring()
def dotplot(
    matrix: Sequence,
    colors: Union[str, Sequence] = "#376B6D",
    xlabels: Optional[Sequence] = None,
    ylabels: Optional[Sequence] = None,
    xlabel_rotation: int = 90,
    ylabel_rotation: int = 0,
    legend_title: Optional[str] = None,
    annotated: bool = False,
    alpha: float = 0.5,
    display: Optional[bool] = None,
    return_plot: bool = False,
    title: Optional[str] = None,
    save: Union[str, Path, None] = None,
):
    """(matplotlib) Dotplot

    Args:
        matrix: Matrix array of values
        colors: The color of every dot
        xlabels: {xlabels}
        ylabels: {ylabels}
        annotated: {annotated}
        xlabel_rotation: {xlabel_rotation}
        ylabel_rotation: {ylabel_rotation}
        legend_title: {legend_title}
        alpha: {alpha}
        title: {title}
        display: {display}
        save: {save}
        return_plot: {return_plot}

    """
    size = np.asarray(matrix)
    (N, M) = size.shape
    X, Y = np.meshgrid(np.arange(M), np.arange(N))
    R = size / size.flatten().max() / 2 * 0.9

    fig, ax = plt.subplots(figsize=(M / 2, N / 2))
    ax.set_aspect("equal")

    circles = [
        plt.Circle((j + 0.5, i + 0.5), radius=r,)
        for r, j, i in zip(R.flat, X.flat, Y.flat)
    ]
    if annotated:
        for v, j, i in zip(size.flat, X.flat, Y.flat):
            if v != 0:
                ax.annotate(v, xy=(j + 0.5, i + 0.5), ha="center", va="center")
    if isinstance(colors, str):
        circ_col = PatchCollection(circles, alpha=alpha, facecolor=colors)
    else:
        color = np.asarray(colors)
        circ_col = PatchCollection(circles, alpha=alpha, facecolor=color.flat)

    ax.add_collection(circ_col)

    ax.set(
        yticks=np.arange(0.5, N + 1),
        xticks=np.arange(0.5, M + 1),
        xticklabels=list(xlabels) + [""],
        yticklabels=list(ylabels) + [""],
    )
    plt.xticks(rotation=xlabel_rotation)
    plt.yticks(rotation=ylabel_rotation)

    for spine in plt.gca().spines.values():
        spine.set_visible(False)

    ax.xaxis.set_ticks_position("none")
    ax.yaxis.set_ticks_position("none")

    legends = [
        Line2D(
            (),
            (),
            color="white",
            marker="o",
            markerfacecolor="black",
            label=f"{int(size.max())}",
            markersize=1.0,
        ),
        Line2D(
            (),
            (),
            color="white",
            marker="o",
            markerfacecolor="black",
            label=f"{int(size.max() / 4 * 3)}",
            markersize=0.75,
        ),
        Line2D(
            (),
            (),
            color="white",
            marker="o",
            markerfacecolor="black",
            label=f"{int(size.max() / 2)}",
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

    if legend_title is None:
        legend_title = "Size"

    legend = ax.legend(
        handles=legends,
        loc="upper left",
        bbox_to_anchor=(1, -0.1, 0, 1),
        title=legend_title,
        markerscale=22,
        labelspacing=1.2,
        frameon=False,
    )
    ax.add_artist(legend)

    if title:
        plt.title(title)

    if save:
        plt.savefig(save, dpi=300, bbox_inches="tight")

    if display is None:
        if CONFIG.WORKING_ENV is None:
            display = False
        else:
            display = True
    if not display:
        plt.close()

    if return_plot:
        return fig, ax
