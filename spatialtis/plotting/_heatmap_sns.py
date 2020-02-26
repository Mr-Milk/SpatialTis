import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

from typing import Union, Sequence, Optional

from .palette import get_linear_colors, get_colors


def heatmap(
        df: pd.DataFrame,
        row_label: Optional[str] = None,
        col_label: Optional[str] = None,
        row_colors: Union[Sequence[str], str, None] = None,
        col_colors: Union[Sequence[str], str, None] = None,
        palette: Union[Sequence[str], str, None] = None,
        colorbar_type: str = 'bar',
        categorical_colorbar_text: Union[Sequence[str], str, None] = None,
        legend_bbox: Optional[Sequence[float]] = None,
        colorbar_bbox: Optional[Sequence[float]] = None,
        **kwargs,
):
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
        row_index = df.index.to_frame(index=False)
        row_index_items = [array for name, array in row_index.iteritems() if name in row_colors]
        for array in row_index_items:
            uni_bars += list(pd.unique(array))

    if col_colors is not None:
        col_index = df.columns.to_frame(index=False)
        col_index_items = [array for name, array in col_index.iteritems() if name in col_colors]
        for array in col_index_items:
            uni_bars += list(pd.unique(array))

    colors = get_colors(len(uni_bars), "Category20")
    colors_bar_mapper = dict(zip(uni_bars, colors))

    if row_colors is not None:
        row_annos = pd.concat([c.map(colors_bar_mapper) for c in row_index_items], axis=1)
        plot_kwargs['row_colors'] = row_annos

    if col_colors is not None:
        col_annos = pd.concat([c.map(colors_bar_mapper) for c in col_index_items], axis=1)
        plot_kwargs['col_colors'] = col_annos

    # handle the colors
    if palette is not None:
        cmap = get_linear_colors(palette)
    else:
        cmap = get_linear_colors("Viridis")

    # handle legend bbox
    default_legend_bbox = (-.25, 0.5)
    default_cbar_bbox = (1.05, .1, .03, .15)
    default_cbar_legend_bbox = (-.25, 0.15)
    if legend_bbox is not None:
        default_legend_bbox = legend_bbox
    if colorbar_bbox is not None:
        default_cbar_bbox = colorbar_bbox
        default_cbar_legend_bbox = colorbar_bbox

    if 'cmap' not in plot_kwargs.keys():
        plot_kwargs['cmap'] = cmap

    # change cbar location
    if 'cbar_pos' not in plot_kwargs.keys():
        plot_kwargs['cbar_pos'] = default_cbar_bbox

    # shrink dendrogram, so that i can place the legend
    if 'dendrogram_ratio' not in plot_kwargs.keys():
        plot_kwargs['dendrogram_ratio'] = 0.035

    # add categorical color bar
    if colorbar_type is 'categorical':
        plot_kwargs['cbar_pos'] = None
        if categorical_colorbar_text is None:
            raise ValueError("'categorical_colorbar_text' should be set if use categorical colorbar")
        else:
            texts = list(categorical_colorbar_text)
            cbar_mapper = zip(texts, [cmap[int(i)] for i in np.linspace(0, len(cmap) - 1, len(texts))])
            cbar_legends = [mpatches.Patch(label=l, color=c) for l, c in cbar_mapper]

    # return a ax_heatmap instance
    h = sns.clustermap(heat_data, **plot_kwargs)

    # add legends for color annotations
    legends = [mpatches.Patch(label=l, color=c) for l, c in colors_bar_mapper.items()]
    ax = h.ax_heatmap

    if len(legends) > 0:
        add_legend = ax.legend(loc="center left",
                               bbox_to_anchor=default_legend_bbox,
                               handlelength=0.8,
                               handles=legends,
                               frameon=False)
        plt.gca().add_artist(add_legend)

    if colorbar_type is 'categorical':
        add_cbar_legend = ax.legend(loc="center left",
                                    bbox_to_anchor=default_cbar_legend_bbox,
                                    handlelength=0.8,
                                    handles=cbar_legends,
                                    frameon=False)

    if row_label is None:
        remove_ytlabel = ax.set_yticklabels("")
        remove_ylabel = ax.set_ylabel("")
        remove_yticks = ax.set_yticks([])

    if col_label is None:
        remove_xtlabel = ax.set_xticklabels("")
        remove_xlabel = ax.set_xlabel("")
        remove_xticks = ax.set_xticks([])
