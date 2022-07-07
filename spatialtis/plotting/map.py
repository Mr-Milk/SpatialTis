from __future__ import annotations

import numpy as np
from anndata import AnnData
from itertools import cycle
from legendkit import CatLegend
from matplotlib import pyplot as plt
from matplotlib.colors import to_hex
from milkviz import point_map, point_map3d, polygon_map
from typing import List, Optional, Tuple, Dict

from spatialtis.abc import AnalysisBase
from spatialtis.utils import doc
from .utils import COLOR_POOL


def _fig_layout(count, ncol):
    if count <= ncol:
        nrow = 1
        ncol = count
    else:
        nrow = count // ncol + (count % ncol > 0)
    return nrow, ncol


def _sep_plot_options(plot_options, hijack="legend_kw"):
    # hijack the legend configuration
    legend_kw = {}
    if hijack in plot_options.keys():
        legend_kw = plot_options[hijack]
        del plot_options[hijack]
    # don't allow use to set ax to prevent overlay
    if 'ax' in plot_options.keys():
        del plot_options['ax']
    return legend_kw, plot_options


def _color_mapper(ab,
                  data,
                  masked_type_name,
                  masked_type_color,
                  types_colors=None,
                  selected_types=None,
                  ):
    # assign cell colors for each cell type
    color_mapper = {}
    legend_color_mapper = {}
    store_key = "cell_colors"
    unique_types = ab.cell_types
    if ab.has_cell_type:
        # if user specific a colormap, we use user one
        # otherwise, create a new one or read from anndata
        if types_colors is not None:
            color_mapper = types_colors
        else:
            if store_key in data.uns_keys():
                # alloc new to prevent modified the stored version
                color_mapper = {**data.uns[store_key]}
            else:
                # Create a global color mapper to ensure that the cell color is the same across ROI
                color_mapper = dict(zip(unique_types, cycle(COLOR_POOL)))
                # alloc new to prevent modified the stored version
                data.uns[store_key] = {**color_mapper}
        if selected_types is not None:
            unique_types = np.unique(selected_types)
            for t in selected_types:
                legend_color_mapper[t] = color_mapper[t]
            masked_type_color = to_hex(masked_type_color, keep_alpha=True)
            color_mapper[masked_type_name] = masked_type_color
            legend_color_mapper[masked_type_name] = masked_type_color
        else:
            legend_color_mapper = {**color_mapper}
    return color_mapper, legend_color_mapper, unique_types


def _masked_cell_type(cell_types, unique_types, selected_types, masked_type_name):
    if (cell_types is not None) & (selected_types is not None):
        cell_mask = np.isin(cell_types, unique_types)
        cell_types[~cell_mask] = masked_type_name
    return cell_types


@doc
def cell_map(
        data: AnnData,
        rois: List[str],
        ncol: int = 5,
        use_shape: bool = False,
        show_neighbors: bool = False,
        selected_types: Optional[List] = None,
        masked_type_name: str = "Other",
        masked_type_color: str = "#d3d3d3",
        figsize: Tuple[float, float] = None,
        wspace: float = 0,
        hspace: float = 0.1,
        types_colors: Dict = None,
        cell_type_key: str = None,
        shape_key: str = None,
        centroid_key: str = None,
        roi_key: str = None,
        **plot_options,
):
    """Visualize cells and neighbors relationship in ROI

    Parameters
    ----------
    data : {adata_plotting}
    rois : list of str
        A list of ROI name that you want to plot.
    ncol : int
        The number of columns in the figure layout.
    use_shape : bool
        Plot cell in polygon only when shape data is available.
    show_neighbors : bool
        Plot the neighbors' relationship.
    selected_types : {selected_types}
    masked_type_name : str, default: 'Other'
        The name of the cell types not in selected_types.
    masked_type_color : color-like, default: '#d3d3d3'
        The color of the cell types not in selected_types.
    figsize : tuple of float
        The size of figure.
    wspace : float, default: 0
        The space between plots vertically.
    hspace : float, default: 0.1
        The space between plots horizontally.
    types_colors : dict
        Change the color for each cell type,
        Key is the type and value is the color.
    cell_type_key : {cell_type_key}
    shape_key : {shape_key}
    centroid_key : {centroid_key}
    roi_key : {roi_key}
    **plot_options:
        Pass to :func:`milkviz.point_map` or :func:`milkviz.point_map3d` or :func:`milkviz.polygon_map`

    """
    ab = AnalysisBase(data,
                      cell_type_key=cell_type_key,
                      shape_key=shape_key,
                      centroid_key=centroid_key,
                      roi_key=roi_key,
                      verbose=False)
    if show_neighbors:
        if ab.dimension == 3:
            raise NotImplementedError("Does not support 3D neighbor map")
        ab.check_neighbors()
    ab.is_rois_name_unique()
    if isinstance(rois, str):
        rois = [rois]

    color_mapper, legend_color_mapper, unique_types = \
        _color_mapper(ab,
                      data,
                      masked_type_name,
                      masked_type_color,
                      types_colors,
                      selected_types
                      )

    roi_count = len(rois)
    nrow, ncol = _fig_layout(roi_count, ncol)
    if figsize is None:
        figsize = (ncol * 4, nrow * 4)
    fig = plt.figure(figsize=figsize)
    legend_kw, plot_options = _sep_plot_options(plot_options)

    ax_index = 1
    axes = []
    if show_neighbors:
        for roi_name, points, cell_types, labels, neighbors in ab.iter_roi(
                fields=['centroid', 'cell_type', 'neighbors'],
                filter_rois=rois,
                disable_pbar=True,
        ):
            cell_types = _masked_cell_type(cell_types, unique_types,
                                           selected_types, masked_type_name)
            # get points
            points = np.array(points)
            x, y = points[:, 0], points[:, 1]
            # get neighbors
            labels = np.asarray(labels)
            nmin = labels.min()
            links = []
            for l, neigh in zip(labels, neighbors):
                for n in neigh:
                    if int(n) > l:
                        links.append((n - nmin, l - nmin))
            ax = fig.add_subplot(nrow, ncol, ax_index)
            point_map(x, y, types=cell_types,
                      links=links,
                      types_colors=color_mapper,
                      ax=ax, legend=False,
                      **plot_options)
            ax.set_title(", ".join([str(i) for i in roi_name]))
            ax_index += 1
            axes.append(ax)
    else:
        for roi_name, points, cell_types, polygons in ab.iter_roi(
                fields=['centroid', 'cell_type', 'shape'],
                filter_rois=rois,
                disable_pbar=True,
        ):
            cell_types = _masked_cell_type(cell_types, unique_types,
                                           selected_types, masked_type_name)
            if use_shape:
                ax = fig.add_subplot(nrow, ncol, ax_index)
                polygon_map(polygons, types=cell_types,
                            types_colors=color_mapper, ax=ax,
                            legend=False, **plot_options)
            else:
                points = np.array(points)
                if ab.dimension == 2:
                    x, y = points[:, 0], points[:, 1]
                    ax = fig.add_subplot(nrow, ncol, ax_index)
                    point_map(x, y, types=cell_types,
                              types_colors=color_mapper, ax=ax,
                              legend=False, **plot_options)
                else:
                    x, y, z = points[:, 0], points[:, 1], points[:, 2]
                    ax = fig.add_subplot(nrow, ncol, ax_index, projection="3d")
                    ax = point_map3d(x, y, z, types=cell_types,
                                     types_colors=color_mapper, ax=ax,
                                     legend=False, **plot_options)
            ax.set_title(", ".join([str(i) for i in roi_name]))
            ax_index += 1
            axes.append(ax)

    if ab.has_cell_type:
        legend_ax = axes[-1]
        labels, colors = zip(*legend_color_mapper.items())
        legend_options = dict(
            title="Cell Type",
            title_align="left",
            bbox_to_anchor=(1.05, 0.5),
            loc="center left",
        )
        legend_options = {**legend_options, **legend_kw}
        CatLegend(colors=colors, labels=labels, handle="circle",
                  ax=legend_ax, **legend_options)
    plt.subplots_adjust(wspace=wspace, hspace=hspace)
    plt.close()
    return fig


@doc
def expression_map(
        data: AnnData,
        rois: List[str],
        markers: List[str],
        use_shape: bool = False,
        x_axis: str = "marker",
        figsize: Tuple = None,
        wspace: float = 0,
        hspace: float = 0.2,
        selected_types: List = None,
        cell_type_key: str = None,
        marker_key: str = None,
        shape_key: str = None,
        centroid_key: str = None,
        roi_key: str = None,
        **plot_options,
):
    """Visualize marker expression in ROI

    Parameters
    ----------
    data : {adata_plotting}
    rois : list of str
        A list of ROI name that you want to plot.
    markers : list of str
        A list of markers name that you want to plot.
    x_axis : {'marker', 'roi'}, default: 'marker'
        What is on the x-axis, the marker or roi.
    use_shape : bool
        Plot cell in polygon only when shape data is available.
    figsize : tuple of float
        The size of figure.
    wspace : float, default: 0
        The space between plots vertically.
    hspace : float, default: 0.1
        The space between plots horizontally.
    selected_types : {selected_types}
    cell_type_key : {cell_type_key}
    marker_key : {marker_key}
    shape_key : {shape_key}
    centroid_key : {centroid_key}
    roi_key : {roi_key}
    **plot_options :
        Pass to :func:`milkviz.point_map` or :func:`milkviz.point_map3d` or :func:`milkviz.polygon_map`

    """
    ab = AnalysisBase(data,
                      cell_type_key=cell_type_key,
                      shape_key=shape_key,
                      centroid_key=centroid_key,
                      roi_key=roi_key,
                      marker_key=marker_key,
                      verbose=False)
    ab.is_rois_name_unique()
    if isinstance(rois, str):
        rois = [rois]
    if isinstance(markers, str):
        markers = [markers]
    unique_types = ab.cell_types
    if ab.has_cell_type & (selected_types is not None):
        unique_types = np.unique(selected_types)

    nrow = len(rois)
    ncol = len(markers)
    ax_indexes = np.arange(1, nrow * ncol + 1)
    if x_axis == "roi":
        nrow, ncol = ncol, nrow
        ax_indexes = ax_indexes.reshape(nrow, ncol).T.flatten()
    if figsize is None:
        figsize = (ncol * 4, nrow * 4)
    fig = plt.figure(figsize=figsize)
    cbar_kw, plot_options = _sep_plot_options(plot_options, "cbar_kw")

    cbar_options = dict(
        orientation="horizontal",
        loc="upper center",
        bbox_to_anchor=(0.5, -0.01)
    )
    cbar_options = {**cbar_options, **cbar_kw}

    ax_indexes_iter = iter(ax_indexes)
    axes = []
    roi_names = []
    for roi_name, points, markers_name, exp, cell_types, polygons in ab.iter_roi(
            fields=['centroid', 'exp', 'cell_type', 'shape'],
            filter_rois=rois,
            disable_pbar=True,
            selected_markers=markers
    ):
        roi_names.append(roi_name)
        cell_mask = None
        if cell_types is not None:
            if selected_types is not None:
                cell_mask = np.isin(cell_types, unique_types)
                exp = exp[:, cell_mask]

        if use_shape:
            if cell_mask is not None:
                polygons = np.asarray(polygons)[cell_mask]
            for varray in exp:
                ax_index = next(ax_indexes_iter)
                ax = fig.add_subplot(nrow, ncol, ax_index)
                polygon_map(polygons, values=varray, ax=ax,
                            cbar_kw=cbar_options, **plot_options)
                axes.append(ax)
        else:
            points = np.array(points)
            if cell_mask is not None:
                points = points[cell_mask]
            if ab.dimension == 2:
                x, y = points[:, 0], points[:, 1]
                for varray in exp:
                    ax_index = next(ax_indexes_iter)
                    ax = fig.add_subplot(nrow, ncol, ax_index)
                    point_map(x, y, values=varray, ax=ax,
                              cbar_kw=cbar_options, **plot_options)
                    axes.append(ax)
            else:
                x, y, z = points[:, 0], points[:, 1], points[:, 2]
                for varray in exp:
                    ax_index = next(ax_indexes_iter)
                    ax = fig.add_subplot(nrow, ncol, ax_index, projection="3d")
                    point_map3d(x, y, z, values=varray, ax=ax,
                                cbar_kw=cbar_options, **plot_options)
                    axes.append(ax)

    # add title

    roi_label = [", ".join([str(i) for i in roi_name]) for roi_name in roi_names]
    label_index = np.arange(0, nrow * ncol).reshape(nrow, ncol)
    x_index = label_index[0]
    y_index = label_index[:, 0]
    x_content = markers
    y_content = roi_label
    if x_axis == "roi":
        axes = np.asarray(axes)[np.argsort(ax_indexes)]
        x_content, y_content = y_content, x_content

    for i, c in zip(x_index, x_content):
        axes[i].set_title(c)
    for i, c in zip(y_index, y_content):
        ax = axes[i]
        if ab.dimension == 2:
            text = getattr(ax, 'text')
        else:
            text = getattr(ax, 'text2D')
        text(-0.01, 0.5, c, transform=ax.transAxes,
             fontdict=dict(rotation=90, va="center", ha="center"))

    plt.subplots_adjust(wspace=wspace, hspace=hspace)
    plt.close()
    return fig
