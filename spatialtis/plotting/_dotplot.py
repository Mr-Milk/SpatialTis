from pathlib import Path
from typing import Optional, Sequence, Union

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.collections import PatchCollection
from matplotlib.lines import Line2D


def dotplot(
    df: pd.DataFrame,
    x: Optional[str] = None,
    y: Optional[str] = None,
    annotated: bool = True,
    xlabel_rotation: int = 90,
    ylabel_rotation: int = 0,
    color: Union[str, Sequence] = "#376B6D",
    alpha: float = 0.5,
    display: bool = True,
    return_plot: bool = False,
    title: Optional[str] = None,
    save: Union[str, Path, None] = None,
):
    xlabels = df.columns.get_level_values(x)
    ylabels = df.index.get_level_values(y)

    (N, M) = df.shape
    X, Y = np.meshgrid(np.arange(M), np.arange(N))
    size = df.to_numpy()
    R = size / size.max() / 2 * 0.9

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
    if isinstance(color, str):
        circ_col = PatchCollection(circles, alpha=alpha, facecolor=color)
    else:
        color = np.array(color)
        circ_col = PatchCollection(circles, alpha=alpha, facecolor=color.flat)

    ax.add_collection(circ_col)

    ax.set(
        yticks=np.arange(0.5, N + 1),
        xticks=np.arange(0.5, M + 1),
        xticklabels=xlabels,
        yticklabels=ylabels,
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

    legend = ax.legend(
        handles=legends,
        loc="upper left",
        bbox_to_anchor=(1, -0.1, 0, 1),
        title="Size",
        markerscale=21,
        labelspacing=1.2,
        frameon=False,
    )
    ax.add_artist(legend)

    if title:
        plt.title(title)

    if save:
        plt.savefig(save, dpi=300, bbox_inches="tight")

    if not display:
        plt.close()

    if return_plot:
        return fig, ax
