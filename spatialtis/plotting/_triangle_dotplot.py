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
    diagnol_arr,
    labels,
    annotate=False,
    xlabel_rotation=90,
    ylabel_rotation=0,
    return_plot=False,
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

    circ_col = PatchCollection(circles, alpha=0.5, facecolor="#2E5C6E")

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

    plt.show()

    if return_plot:
        return fig, ax
