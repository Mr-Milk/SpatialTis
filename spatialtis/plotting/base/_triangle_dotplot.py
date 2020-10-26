from pathlib import Path
from typing import Optional, Sequence, Union

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import PatchCollection
from matplotlib.lines import Line2D

from spatialtis.config import CONFIG
from spatialtis.utils import reuse_docstring

# the relationship between label and arr should look like this

"""
    A B C D
A   1 2 3 4
B   1 2 3
C   1 2
D   1
"""


@reuse_docstring()
def tri_dotplot(
    diagonal_arr: Sequence,
    colors: Union[str, Sequence] = "#376B6D",
    labels: Optional[Sequence] = None,
    xlabel_rotation: int = 90,
    ylabel_rotation: int = 0,
    legend_title: Optional[str] = None,
    annotate: bool = False,
    alpha: float = 0.5,
    display: Optional[bool] = None,
    title: Optional[str] = None,
    save: Union[str, Path, None] = None,
    return_plot: bool = False,
):
    """(matplotlib) Triangular dotplot

    This input diagnol array should look like this

    .. code-block:: none

               A B C D
            A  1 2 3 4
            B  1 2 3
            C  1 2
            D  1

    Args:
        diagonal_arr: The diagonal array
        labels: Labels will be marked in x and y axis
        colors: Colors of the dot
        xlabel_rotation: {xlabel_rotation}
        ylabel_rotation: {ylabel_rotation}
        legend_title: {legend_title}
        annotate: {annotate}
        alpha: {alpha}
        display: {display}
        title: {title}
        save: {save}
        return_plot: {return_plot}

    """

    M = N = len(labels)
    flat_b = []
    for arr in diagonal_arr:
        for i in arr:
            flat_b.append(i)
    cm = max(flat_b)

    fig, ax = plt.subplots(figsize=(N / 2, M / 2))
    ax.set_aspect("equal")

    circles = list()

    for i, arr in enumerate(diagonal_arr):
        for t, v in enumerate(arr[::-1]):
            if v != 0:
                circles.append(plt.Circle((i + 0.5, t + 0.5), radius=v / cm / 2 * 0.9))
                if annotate:
                    ax.annotate(v, xy=(i + 0.5, t + 0.5), ha="center", va="center")
    if isinstance(colors, str):
        circ_col = PatchCollection(circles, alpha=alpha, facecolor=colors)
    else:
        color = np.asarray(colors)
        circ_col = PatchCollection(circles, alpha=alpha, facecolor=color.flat)

    ax.add_collection(circ_col)

    ax.set(
        yticks=np.arange(0.5, M + 1),
        xticks=np.arange(0.5, N + 1),
        xticklabels=list(labels) + [""],
        yticklabels=list(labels)[::-1] + [""],
    )

    legends = [
        Line2D(
            (),
            (),
            color="white",
            marker="o",
            markerfacecolor="black",
            label=f"{int(cm)}",
            markersize=1.0,
        ),
        Line2D(
            (),
            (),
            color="white",
            marker="o",
            markerfacecolor="black",
            label=f"{int(cm / 4 * 3)}",
            markersize=0.75,
        ),
        Line2D(
            (),
            (),
            color="white",
            marker="o",
            markerfacecolor="black",
            label=f"{int(cm / 2)}",
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

    plt.xticks(rotation=xlabel_rotation)
    plt.yticks(rotation=ylabel_rotation)

    for spine in plt.gca().spines.values():
        spine.set_visible(False)

    ax.xaxis.set_ticks_position("none")
    ax.yaxis.set_ticks_position("none")

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
