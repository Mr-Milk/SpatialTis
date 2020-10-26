from pathlib import Path
from typing import Optional, Sequence, Union

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from spatialtis.config import CONFIG
from spatialtis.utils import reuse_docstring

from .palette import get_colors, get_linear_colors


@reuse_docstring()
def heatmap(
    df: pd.DataFrame,
    row_label: Optional[str] = None,
    col_label: Optional[str] = None,
    row_colors: Union[Sequence[str], str, None] = None,
    col_colors: Union[Sequence[str], str, None] = None,
    palette: Union[Sequence[str], str, None] = None,
    colorbar_type: str = "bar",
    categorical_colorbar_text: Union[Sequence[str], str, None] = None,
    row_colors_legend_bbox: Optional[Sequence[float]] = None,
    col_colors_legend_bbox: Optional[Sequence[float]] = None,
    colorbar_bbox: Optional[Sequence[float]] = None,
    save: Union[Path, str, None] = None,
    display: Optional[bool] = None,
    return_plot: bool = False,
    **kwargs,
):
    """(matplotlib) A higher wrapper for seaborn's clustermap

        The order of cateforical_colorbar_text match to the order of palette

    Args:
        df: Input data, all data should store in data part, annotations should store in MultiIndex
        row_label: Which level to plot row text
        col_label: Which level to plot col text
        row_colors: Which level to plot row colors
        col_colors: Which level to plot col colors
        palette: {palette}
        colorbar_type: Options are 'continuous' and 'categorical'
        categorical_colorbar_text: If 'categorical', the text need to be provided
        row_colors_legend_bbox: Adjust the locations of row colors legend
        col_colors_legend_bbox: Adjust the locations of col colors legend
        colorbar_bbox: Adjust the locations of colorbar
        display: {display}
        save: {save}
        return_plot: {return_plot}

    """
    try:
        plot_kwargs = dict(**kwargs)
    except NameError:
        plot_kwargs = dict()
    heat_data = pd.DataFrame(df.to_numpy())
    # ==============handle axis test====================
    if row_label is not None:
        row_index = df.index.to_frame(index=False)
        heat_data.index = row_index[row_label]

    if col_label is not None:
        col_index = df.columns.to_frame(index=False)
        heat_data.columns = col_index[col_label]

    # ==============handle color labels====================
    uni_bars = []
    if row_colors is not None:
        row_colors_legend = []
        row_index = df.index.to_frame(index=False)
        row_index_items = [
            array for name, array in row_index.iteritems() if name in row_colors
        ]
        for array in row_index_items:
            row_colors_legend += list(pd.unique(array))
        row_colors_legend = list(pd.unique(row_colors_legend))
        uni_bars += row_colors_legend

    if col_colors is not None:
        col_colors_legend = []
        col_index = df.columns.to_frame(index=False)
        col_index_items = [
            array for name, array in col_index.iteritems() if name in col_colors
        ]
        for array in col_index_items:
            col_colors_legend += list(pd.unique(array))
        col_colors_legend = list(pd.unique(col_colors_legend))
        uni_bars += col_colors_legend

    colors = get_colors(len(uni_bars), ["Category20", "Category20b"])
    colors_bar_mapper = dict(zip(uni_bars, colors))

    if row_colors is not None:
        # print(colors_bar_mapper)
        # print(row_index_items)
        row_annos = pd.concat(
            [c.map(colors_bar_mapper) for c in row_index_items], axis=1
        )
        plot_kwargs["row_colors"] = row_annos

    if col_colors is not None:
        col_annos = pd.concat(
            [c.map(colors_bar_mapper) for c in col_index_items], axis=1
        )
        plot_kwargs["col_colors"] = col_annos

    # handle the colors
    if palette is not None:
        cmap = get_linear_colors(palette)
    else:
        cmap = get_linear_colors(["Viridis"])

    # handle legend bbox

    default_cbar_bbox = (1.05, 0.1, 0.03, 0.15)

    default_col_legend_bbox = (-0.25, 0.85)
    default_row_legend_bbox = (-0.25, 0.5)
    default_cbar_legend_bbox = (-0.25, 0.15)
    if row_colors_legend_bbox is not None:
        default_row_legend_bbox = row_colors_legend_bbox
    if col_colors_legend_bbox is not None:
        default_col_legend_bbox = col_colors_legend_bbox
    if colorbar_bbox is not None:
        default_cbar_bbox = colorbar_bbox
        default_cbar_legend_bbox = colorbar_bbox

    if "cmap" not in plot_kwargs.keys():
        plot_kwargs["cmap"] = cmap

    # change cbar location
    if "cbar_pos" not in plot_kwargs.keys():
        plot_kwargs["cbar_pos"] = default_cbar_bbox

    # shrink dendrogram, so that i can place the legend
    if "dendrogram_ratio" not in plot_kwargs.keys():
        plot_kwargs["dendrogram_ratio"] = 0.035

    # add categorical color bar
    if colorbar_type == "categorical":
        plot_kwargs["cbar_pos"] = None
        if categorical_colorbar_text is None:
            raise ValueError(
                "'categorical_colorbar_text' should be set if use categorical colorbar"
            )
        else:
            texts = categorical_colorbar_text
            cbar_mapper = zip(
                texts, [cmap[int(i)] for i in np.linspace(0, len(cmap) - 1, len(texts))]
            )
            cbar_legends = [mpatches.Patch(label=l, color=c) for l, c in cbar_mapper]

    # return a ax_heatmap instance
    h = sns.clustermap(heat_data, **plot_kwargs)
    ax = h.ax_heatmap

    # add legends for color annotations

    if row_colors is not None:
        row_legends = [
            mpatches.Patch(label=l, color=colors_bar_mapper[l])
            for l in row_colors_legend
        ]
        add_row_legend = ax.legend(
            loc="center left",
            bbox_to_anchor=default_row_legend_bbox,
            handlelength=0.8,
            handles=row_legends,
            frameon=False,
        )
        ax.add_artist(add_row_legend)

    if col_colors is not None:
        col_legends = [
            mpatches.Patch(label=l, color=colors_bar_mapper[l])
            for l in col_colors_legend
        ]
        add_col_legend = ax.legend(
            loc="center left",
            bbox_to_anchor=default_col_legend_bbox,
            handlelength=0.8,
            handles=col_legends,
            frameon=False,
        )
        ax.add_artist(add_col_legend)

    if colorbar_type == "categorical":
        add_cbar_legend = ax.legend(
            loc="center left",
            bbox_to_anchor=default_cbar_legend_bbox,
            handlelength=0.8,
            handles=cbar_legends,
            frameon=False,
        )

    if row_label is None:
        remove_ytlabel = ax.set_yticklabels("")
        remove_ylabel = ax.set_ylabel("")
        remove_yticks = ax.set_yticks([])

    if col_label is None:
        remove_xtlabel = ax.set_xticklabels("")
        remove_xlabel = ax.set_xlabel("")
        remove_xticks = ax.set_xticks([])

    if save:
        h.savefig(save, dpi=300)

    if display is None:
        if CONFIG.WORKING_ENV is None:
            display = False
        else:
            display = True
    if not display:
        plt.close()

    if return_plot:
        return h
