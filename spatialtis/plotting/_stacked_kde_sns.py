from typing import Sequence, Union

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from .palette import get_colors, get_linear_colors


def stacked_kde(
    df: pd.DataFrame,
    col: str,
    row: str,
    palette: Union[Sequence[str], str, None] = None,
    display: bool = True,
    return_plot: bool = False,
):
    new_df = pd.DataFrame(df.stack(), columns=["value"]).reset_index()

    # handle colors
    default_palette = get_colors(9, ["Set1"])
    if palette is not None:
        default_palette = get_linear_colors(palette)

    colormapper = dict(zip(pd.unique(new_df[col]), default_palette))

    def kde(x, t, **kwargs):
        name = np.unique(t)[0]
        c = colormapper[name]
        mean = np.mean(x)
        ax = sns.kdeplot(x, color=c, shade=True)
        ylim = ax.get_ylim()
        plt.vlines(mean, ylim[0], ylim[1], color=c)

    g = sns.FacetGrid(new_df, col=col, row=row, margin_titles=True)
    g.map(kde, "value", col)

    for ax in g.axes.flat:
        plt.setp(ax.texts, text="")
        ax.set_xlabel("")
        ax.set_ylabel("")

    g.set_titles(col_template="{col_name}", row_template="{row_name}")

    if not display:
        plt.close()

    if return_plot:
        return g
