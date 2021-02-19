from typing import Union

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.collections import PatchCollection
from matplotlib.lines import Line2D

from spatialtis.plotting.abc import MatplotlibMixin
from spatialtis.typing import Array
from spatialtis.utils import doc


@doc
class dot_plot(MatplotlibMixin):
    """Dot plot, Matplotlib

    Args:
        data: Matrix array of values
        colors: The color of every dot
        annotated: To show value number on the dot
        alpha: Alpha value for opacity
        **plot_options: {plot_options}

    """

    def __init__(
        self,
        data: pd.DataFrame,
        colors: Union[str, Array] = "#376B6D",
        annotated: bool = False,
        alpha: float = 0.5,
        **plot_options,
    ):
        super().__init__(**plot_options)
        self.data = data
        data_max = data.max().max()
        (N, M) = data.shape
        X, Y = np.meshgrid(np.arange(M), np.arange(N))
        R = data / data_max / 2 * 0.9
        if self.size is None:
            self.size = (
                M / 2,
                N / 2,
            )
        self.fig, self.ax = plt.subplots(figsize=self.size)
        self.ax.set_aspect("equal")

        circles = [
            plt.Circle((j + 0.5, i + 0.5), radius=r,)
            for r, j, i in zip(R.to_numpy().flat, X.flat, Y.flat)
        ]
        if annotated:
            for v, j, i in zip(data.to_numpy().flat, X.flat, Y.flat):
                if (v != 0) & (not np.isnan(v)):
                    self.ax.annotate(v, xy=(j + 0.5, i + 0.5), ha="center", va="center")
        if isinstance(colors, str):
            circ_col = PatchCollection(circles, alpha=alpha, facecolor=colors)
        else:
            color = np.asarray(colors)
            circ_col = PatchCollection(circles, alpha=alpha, facecolor=color.flat)

        self.ax.add_collection(circ_col)

        self.ax.set(
            yticks=np.arange(0.5, N + 1),
            xticks=np.arange(0.5, M + 1),
            xticklabels=list(data.columns) + [""],
            yticklabels=list(data.index) + [""],
        )
        plt.xticks(rotation=self.xtickslabel_rotation)
        plt.yticks(rotation=self.ytickslabel_rotation)

        # for spine in plt.gca().spines.values():
        #     spine.set_visible(False)
        #
        self.ax.xaxis.set_ticks_position("none")
        self.ax.yaxis.set_ticks_position("none")
        self.ax.grid(False)
        legends = [
            Line2D(
                (),
                (),
                color="white",
                marker="o",
                markerfacecolor="black",
                label=f"{int(data_max * ratio)}",
                markersize=ratio,
            )
            for ratio in [1.0, 0.75, 0.5, 0.25]
        ]

        if self.legend_title is None:
            self.legend_title = "Size"

        legend = self.ax.legend(
            handles=legends,
            loc="upper left",
            bbox_to_anchor=(1, -0.1, 0, 1),
            title=self.legend_title,
            markerscale=22,
            labelspacing=1.2,
            frameon=False,
        )
        self.ax.add_artist(legend)
        self.set_up()
