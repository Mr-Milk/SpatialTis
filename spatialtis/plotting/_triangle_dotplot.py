from pathlib import Path
from typing import Optional, Sequence, Union

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import PatchCollection

# the relationship between label and arr should look like this

"""
    A B C D
A   1 2 3 4
B   1 2 3
C   1 2
D   1
"""


def tri_dotplot(
    diagnol_arr: np.array,
    labels: Sequence,
    annotate: bool = False,
    xlabel_rotation: int = 90,
    ylabel_rotation: int = 0,
    color: Union[str, Sequence] = "#376B6D",
    return_plot: bool = False,
    display: bool = True,
    title: Optional[str] = None,
    save: Union[str, Path, None] = None,
):
    all_max = list()
    for arr in diagnol_arr:
        all_max.append(max(arr))

    M = N = len(labels)
    cm = max(all_max)

    fig, ax = plt.subplots(figsize=(N / 2, M / 2))

    circles = list()

    for i, arr in enumerate(diagnol_arr):
        for t, v in enumerate(arr[::-1]):
            if v != 0:
                circles.append(plt.Circle((i + 0.5, t + 0.5), radius=v / cm / 2 * 0.9))
                if annotate:
                    ax.annotate(v, xy=(i + 0.5, t + 0.5), ha="center", va="center")

    circ_col = PatchCollection(circles, alpha=0.5, facecolor=color)

    ax.add_collection(circ_col)

    ax.set(
        yticks=np.arange(0.5, M + 1),
        xticks=np.arange(0.5, N + 1),
        xticklabels=labels,
        yticklabels=labels[::-1],
    )
    plt.xticks(rotation=xlabel_rotation)
    plt.yticks(rotation=ylabel_rotation)

    for spine in plt.gca().spines.values():
        spine.set_visible(False)

    ax.xaxis.set_ticks_position("none")
    ax.yaxis.set_ticks_position("none")

    if title:
        plt.title(title)

    if save:
        fig.savefig(save, dpi=300)

    if not display:
        plt.close()

    if return_plot:
        return fig, ax
