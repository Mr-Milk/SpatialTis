from typing import Dict, Optional, Sequence, Union

import matplotlib.patches as mpatches
import numpy as np
import pandas as pd
import seaborn as sns

from spatialtis.utils import doc

from ..abc import MatplotlibMixin
from .palette import get_colors, get_linear_colors


@doc
class heatmap(MatplotlibMixin):
    """Heatmap, Matplotlib

    The order of cateforical_colorbar_text match to the order of palette

    Args:
        df: Input data, all data should store in data part, annotations should store in MultiIndex
        row_label: Which level to plot row text
        col_label: Which level to plot col text
        row_colors: Which level to plot row colors
        col_colors: Which level to plot col colors
        categorical_colorbar: Colorbar in categorical way, the text need to be provided
        clustermap_kwargs: Pass to `seaborn.clustermap <https://seaborn.pydata.org/generated/seaborn.clustermap.html>`_
        **plot_options: {plot_options}

    """

    def __init__(
        self,
        df: pd.DataFrame,
        row_label: Optional[str] = None,
        col_label: Optional[str] = None,
        row_colors: Union[Sequence[str], str, None] = None,
        col_colors: Union[Sequence[str], str, None] = None,
        categorical_colorbar: Optional[Sequence[str]] = None,
        clustermap_kwargs: Optional[Dict] = None,
        **plot_options,
    ):
        super().__init__(**plot_options)
        if clustermap_kwargs is None:
            clustermap_kwargs = {}
        raw = df.to_numpy()
        raw_min = raw.min()
        raw_max = raw.max()
        plot_kwargs = {
            "method": "ward",
            "cbar_pos": (1, 0.1, 0.05, 0.15),
            "cbar_kws": dict(ticks=np.linspace(raw_min, raw_max, num=4)),
            **clustermap_kwargs,
        }
        heat_data = pd.DataFrame(df.to_numpy())
        self.row_index = df.index.to_frame(index=False)
        self.col_index = df.columns.to_frame(index=False)
        # ==============handle axis test====================
        if row_label is not None:
            heat_data.index = self.row_index[row_label]

        if col_label is not None:
            heat_data.columns = self.col_index[col_label]

        # ==============handle color labels====================
        colors_pool = get_colors(60, ["Category20", "Category20b", "Set1", "Paired"])
        if row_colors is not None:
            row_colors_labels = pd.unique(
                self.row_index[row_colors].to_numpy().T.flatten()
            )
            row_colors_mapper = dict(
                zip(row_colors_labels, colors_pool[0 : len(row_colors_labels)])
            )
            colors_pool = colors_pool[len(row_colors_labels) : :]
            row_annos = self.row_index[row_colors].replace(row_colors_mapper)
            plot_kwargs["row_colors"] = row_annos

        if col_colors is not None:
            col_colors_labels = pd.unique(
                self.col_index[col_colors].to_numpy().T.flatten()
            )
            col_colors_mapper = dict(
                zip(col_colors_labels, colors_pool[0 : len(col_colors_labels)])
            )
            col_annos = self.col_index[col_colors].replace(col_colors_mapper)
            plot_kwargs["col_colors"] = col_annos

        # handle the colors
        default_palette = ["Viridis"]
        if self.palette is None:
            cmap = get_linear_colors(default_palette)
        else:
            cmap = get_linear_colors(self.palette)
        if "cmap" not in plot_kwargs.keys():
            plot_kwargs["cmap"] = cmap

        # shrink dendrogram, so that i can place the legend
        if "dendrogram_ratio" not in plot_kwargs.keys():
            plot_kwargs["dendrogram_ratio"] = 0.035

        # add categorical color bar
        if categorical_colorbar is not None:
            plot_kwargs["cbar_pos"] = None
            cbar_mapper = zip(
                categorical_colorbar,
                [
                    cmap[int(i)]
                    for i in np.linspace(0, len(cmap) - 1, len(categorical_colorbar))
                ],
            )
            cbar_legends = [mpatches.Patch(label=l, color=c) for l, c in cbar_mapper]

        # return a ax_heatmap instance
        self.fig = sns.clustermap(heat_data, **plot_kwargs)

        # add legends for color annotations
        if row_colors is not None:
            row_legends = [
                mpatches.Patch(label=l, color=row_colors_mapper[l])
                for l in row_colors_mapper
            ]
            self.fig.ax_row_dendrogram.legend(
                loc="lower left",
                borderaxespad=0,
                bbox_to_anchor=(-3, 0.5),
                handlelength=0.8,
                handles=row_legends,
                frameon=False,
            )

        if col_colors is not None:
            col_legends = [
                mpatches.Patch(label=l, color=col_colors_mapper[l])
                for l in col_colors_labels
            ]
            self.fig.ax_col_dendrogram.legend(
                loc="lower center",
                borderaxespad=0,
                bbox_to_anchor=(0.5, 1.05),
                ncol=int(len(col_legends) / 2),
                handlelength=0.8,
                handles=col_legends,
                frameon=False,
            )

        if categorical_colorbar is not None:
            self.fig.ax_heatmap.legend(
                loc="lower left",
                borderaxespad=0,
                bbox_to_anchor=(1.1, 0.05),
                handlelength=0.8,
                handles=cbar_legends,
                frameon=False,
            )

        if row_label is None:
            self.fig.ax_heatmap.set_yticklabels("")
            self.fig.ax_heatmap.set_ylabel("")
            self.fig.ax_heatmap.set_yticks([])

        if col_label is None:
            self.fig.ax_heatmap.set_xticklabels("")
            self.fig.ax_heatmap.set_xlabel("")
            self.fig.ax_heatmap.set_xticks([])

        self.set_up()
