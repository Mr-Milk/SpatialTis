import numpy as np
import pandas as pd

from pathlib import Path
import shutil

from svgutils.compose import Figure, SVG

from scipy.spatial.distance import pdist
from scipy.cluster.hierarchy import linkage, dendrogram

from bokeh.models import ColumnDataSource, LinearColorMapper, ColorBar, Legend, LegendItem
from bokeh.plotting import figure, output_file, show
from bokeh.io import output_notebook, export_svgs
from bokeh.layouts import gridplot

from spatialtis.plotting import get_linear_colors, get_colors
from typing import Optional, Sequence, Union


def dendrogram_coord(matrix, metric='euclidean', method='ward'):
    X = pdist(matrix, metric=metric)
    Z = linkage(X, method)
    results = dendrogram(Z, no_plot=True)
    labels = list(map(int, results['ivl']))
    icd, dcd = results['icoord'], results['dcoord']
    return icd, dcd, labels


def heatmap(
        df: pd.DataFrame,
        z_score: Optional[str] = None,
        row_cluster: bool = True,
        col_cluster: bool = True,
        row_label: Optional[str] = None,
        col_label: Optional[str] = None,
        row_colors: Union[Sequence[str], str, None] = None,
        col_colors: Union[Sequence[str], str, None] = None,

        row_label_show: bool = True,
        col_label_show: bool = True,
        label_standoff: float = 1.0,
        colors_standoff: float = 0,
        row_colors_show: bool = True,
        col_colors_show: bool = True,
        row_colors_legend: bool = True,
        col_colors_legend: bool = True,
        row_colors_legend_title: Optional[str] = None,
        col_colors_legend_title: Optional[str] = None,
        row_dendrogram: bool = True,
        col_dendrogram: bool = True,
        color_bar: bool = True,

        dendrogram_ratio: float = 0.1,
        colors_size: float = 10.0,
        method: str = 'ward',
        metric: str = 'euclidean',

        display: bool = True,
        title: Optional[str] = None,
        size: Optional[Sequence[int]] = None,
        save_svg: Optional[str] = None,
        save_html: Optional[str] = None,
        palette: Union[Sequence[str], str, None] = None,
        return_plot: bool = False,
):
    df = df.loc[::-1]

    col_depth = df.columns.nlevels
    index_depth = df.index.nlevels

    col_names = df.columns.names
    index_names = df.index.names

    # assign annotations to row/col label and color-rect
    row_all, col_all = [], []
    if row_label is not None:
        row_all.append(row_label)
    else:
        row_label_show = False

    if col_label is not None:
        col_all.append(col_label)
    else:
        col_label_show = False

    if isinstance(row_colors, str):
        if row_colors not in row_all:
            row_all.append(row_colors)
        else:
            raise ValueError(f"Replicated key in row_label and row_colors")
        row_colors = [row_colors]
    elif isinstance(row_colors, Sequence):
        row_all += row_colors
    else:
        row_colors_show = False
        row_colors_legend = False

    if isinstance(col_colors, str):
        if col_colors not in col_all:
            col_all.append(col_colors)
        else:
            raise ValueError(f"Replicated key in col_label and col_colors")
        col_colors = [col_colors]
    elif isinstance(col_colors, Sequence):
        col_all += col_colors
    else:
        col_colors_show = False
        col_colors_legend = False

    # check if match to the input df
    for r in row_all:
        if r not in index_names:
            raise ValueError(f"{r} not found in rows' names")

    for c in col_all:
        if c not in col_names:
            raise ValueError(f"{c} not found in columns' names")

    # initiate the plot instance
    TOOLS = ""

    figure_config = dict(
        title=title,
        tools=TOOLS,
        toolbar_location=None,
        y_axis_location="right")
    if size is not None:
        figure_config['plot_width'] = size[0]
        figure_config['plot_height'] = size[1]

    p = figure(**figure_config)
    lg = figure(plot_height=int(p.plot_height), plot_width=200, tools="")

    # calculate row_dendrogram
    if row_cluster:
        r_icd, r_dcd, r_labels = dendrogram_coord(df.values, metric=metric, method=method)
        df = df.iloc[r_labels]

    # calculate col_dendrogram
    if col_cluster:
        c_icd, c_dcd, c_labels = dendrogram_coord(df.T.values, metric=metric, method=method)
        df = df.iloc[:, c_labels]

    # if z_score is set
    if z_score == "column":
        df = (df - df.mean(axis=0)) / df.std(axis=0)
    elif z_score == "row":
        df = df.sub(df.mean(axis=1), axis=0).div(df.std(axis=1), axis=0)

    default_palette = ['Cividis', 'Gray', 'Inferno', 'Magma', 'Viridis']
    if isinstance(palette, str):
        if palette in default_palette:
            colors = get_linear_colors(256, palette)
        else:
            raise ValueError(f"Not a build-in palette, Options are 'Cividis', 'Gray', 'Inferno', 'Magma', 'Viridis'")
    elif palette is None:
        colors = get_linear_colors(50, [default_palette[-1]])
    else:
        colors = palette

    nd = df.to_numpy()
    try:
        mapper = LinearColorMapper(palette=colors, low=np.min(nd), high=np.max(nd))
    except ValueError:
        raise ValueError(f"Plotting heatmap only accept built-in linear palette or Sequence of user-defined colors")

    x_stepper = np.arange(0.5, len(df.columns) + 0.5)
    y_stepper = np.arange(0.5, len(df.index) + 0.5)
    xsmax = np.max(x_stepper)
    xsmin = np.min(x_stepper)
    ysmax = np.max(y_stepper)
    ysmin = np.min(y_stepper)

    col_coord = np.array([x_stepper for i in range(0, len(df.index))]).ravel()
    index_coord = np.array([y_stepper for i in range(0, len(df.columns))]).T.ravel()

    data = pd.DataFrame(dict(
        col_coord=col_coord,
        index_coord=index_coord,
        values=nd.ravel(),
    ))

    source = ColumnDataSource(data)

    # plot the heatmap main part
    p.rect(x='col_coord', y='index_coord', width=1, height=1,
           source=source,
           fill_color={'field': 'values', 'transform': mapper},
           line_color=None)

    # we need to calculate the relative pixel level base on the canvas size
    ph = p.plot_height
    pw = p.plot_width
    min_side = np.min([ph, pw])
    y_len = len(y_stepper)
    x_len = len(x_stepper)
    x_scale = x_len / pw
    y_scale = y_len / ph
    # for row_block_width and col_block_height
    row_block_width = colors_size * x_scale if row_colors_show else 0
    col_block_height = colors_size * y_scale if col_colors_show else 0
    row_colors_standoff = colors_standoff * x_scale if row_colors_show else 0
    col_colors_standoff = colors_standoff * y_scale if col_colors_show else 0
    row_label_standoff = label_standoff * x_scale
    col_label_standoff = label_standoff * y_scale

    row_dendro_offset = row_colors_standoff + row_block_width * len(row_colors) if row_colors_show else 0
    col_dendro_offset = col_colors_standoff + col_block_height * len(col_colors) if col_colors_show else 0

    # plot row-dendrogram
    if row_cluster & row_dendrogram:
        row_icd = (r_icd - np.min(r_icd)) / (np.max(r_icd) - np.min(r_icd)) * (ysmax - ysmin) + 0.5
        row_dcd = (r_dcd - np.min(r_dcd)) / (
                    np.max(r_dcd) - np.min(r_dcd)) * dendrogram_ratio * x_len + row_dendro_offset

        p.multi_line(list(-row_dcd), list(row_icd), line_alpha=0.8, line_color='black')

    # plot col-dendrogram
    if col_cluster & col_dendrogram:
        col_icd = (c_icd - np.min(c_icd)) / (np.max(c_icd) - np.min(c_icd)) * (xsmax - xsmin) + 0.5
        col_dcd = (c_dcd - np.min(c_dcd)) / (
                    np.max(c_dcd) - np.min(c_dcd)) * dendrogram_ratio * y_len + ysmax + 0.5 + col_dendro_offset

        p.multi_line(list(col_icd), list(col_dcd), line_alpha=0.8, line_color='black')

    # plot row-colors
    if row_colors_show:
        if index_depth == 1:
            row_blocks = np.array(df.index)
        else:
            df_index = pd.DataFrame([list(i) for i in df.index], columns=df.index.names)
            row_blocks = df_index[row_colors].to_numpy()

        unigp = np.array([np.unique(b) for b in row_blocks.T[::-1]]).ravel()
        row_blocks_mapper = dict(zip(unigp, get_colors(len(unigp), 'Category20')))

        def add_row_blocks(blocks, step):
            plot_blocks = ColumnDataSource(dict(
                x=np.array([-row_block_width / 2 - row_block_width * step for i in
                            range(0, len(blocks))]) - row_colors_standoff,
                y=y_stepper,
                colors=[row_blocks_mapper[n] for n in blocks],
                name=[n for n in blocks]
            ))
            p.rect(x='x', y='y', width=row_block_width, height=1, source=plot_blocks, fill_color='colors',
                   line_color=None)

        if index_depth > 1:
            for step, b in enumerate(row_blocks.T[::-1]):
                add_row_blocks(b, step)
        else:
            add_row_blocks(row_blocks, 0)

        # add row_colors_legends
        if row_colors_legend:
            row_colors_legends = []
            # create invisible glyph rect, so we can add to legends
            for name, color in row_blocks_mapper.items():
                a = lg.rect(1, 1, 1, 1, fill_color=color, line_color=None)
                row_colors_legends.append(LegendItem(label=name, renderers=[a]))

            row_legend = Legend(title=row_colors_legend_title, items=row_colors_legends, location='center_left',
                                label_standoff=0,
                                border_line_color=None)
            lg.add_layout(row_legend, 'center')

    # plot col-colors
    if col_colors_show:
        if col_depth == 1:
            col_blocks = np.array(df.columns)
        else:
            df_col = pd.DataFrame([list(i) for i in df.columns], columns=df.columns.names)
            col_blocks = df_col[col_colors].to_numpy()

        unigp = np.array([np.unique(b) for b in col_blocks.T[::-1]]).ravel()
        col_blocks_mapper = dict(zip(unigp, get_colors(len(unigp), 'Spectral')))
        col_colors_legends = []

        def add_col_blocks(blocks, step):
            plot_blocks = ColumnDataSource(dict(
                x=x_stepper,
                y=np.array([ysmax + 0.5 + col_block_height / 2 + col_block_height * step for i in
                            range(0, len(blocks))]) + col_colors_standoff,
                colors=[col_blocks_mapper[n] for n in blocks],
                name=[n for n in blocks]
            ))
            p.rect(x='x', y='y', height=col_block_height, width=1, source=plot_blocks, fill_color='colors',
                   line_color=None)

        if col_depth > 1:
            for step, b in enumerate(col_blocks.T[::-1]):
                add_col_blocks(b, step)
        else:
            add_col_blocks(col_blocks, 0)

        # add col_colors_legends
        if col_colors_legend:
            col_colors_legends = []
            # create invisible glyph rect, so we can add to legends
            for name, color in col_blocks_mapper.items():
                a = lg.rect(1, 1, 1, 1, fill_color=color, line_color=None)
                col_colors_legends.append(LegendItem(label=name, renderers=[a]))

            col_legend = Legend(title=col_colors_legend_title, items=col_colors_legends, location="top_left",
                                label_standoff=0,
                                border_line_color=None)
            lg.add_layout(col_legend, 'center')

    # plot color bar
    if color_bar:
        color_bar = ColorBar(color_mapper=mapper, major_label_text_font_size="8pt",
                             # ticker=BasicTicker(desired_num_ticks=5),
                             # bar_line_color=None,
                             # formatter=PrintfTickFormatter(format="%d%%"),
                             height=100,
                             label_standoff=6,
                             border_line_color=None, location=(0, 0))
        lg.add_layout(color_bar, 'center')

    # draw a white rect as legend background
    lg.rect(1, 1, 1, 1, fill_color="white", line_color="white")

    # plot row-label
    p.axis.major_label_text_color = None
    if row_label_show:
        if index_depth == 1:
            row_text = list(df.index)
        else:
            df_index = pd.DataFrame([list(i) for i in df.index], columns=df.index.names)
            row_text = df_index[row_label].to_numpy()

        p.yaxis.major_label_text_color = 'black'
        p.yaxis.fixed_location = xsmax + 0.5
        p.yaxis.major_label_standoff = - 5
        p.yaxis.ticker = y_stepper
        p.yaxis.major_label_overrides = {k: v for k, v in zip(y_stepper, row_text)}

    if col_label_show:
        if col_depth == 1:
            col_text = list(df.columns)
        else:
            df_col = pd.DataFrame([list(i) for i in df.columns], columns=df.columns.names)
            col_text = df_col[col_label].to_numpy()

        p.xaxis.fixed_location = ysmin - 0.5
        p.yaxis.ticker = x_stepper
        p.yaxis.major_label_overrides = {k: v for k, v in zip(x_stepper, col_text)}

    # doing some adjustments
    p.grid.grid_line_color = None
    p.axis.axis_line_color = None
    p.axis.major_tick_line_color = None
    p.axis.minor_tick_line_color = None
    p.axis.major_label_text_align = "center"
    p.axis.major_label_text_baseline = "middle"
    p.hover.point_policy = "follow_mouse"
    # p.outline_line_color = None

    lg.grid.grid_line_color = None
    lg.axis.axis_line_color = None
    lg.axis.major_tick_line_color = None
    lg.axis.minor_tick_line_color = None
    lg.axis.major_label_text_color = None
    # lg.outline_line_color = None

    if row_colors_legend | col_colors_legend:
        lg.legend.label_text_font_size = '10pt'
        lg.legend.glyph_width = 15
        lg.legend.glyph_height = 15
        lg.legend.label_text_baseline = "middle"
        lg.legend.label_standoff = 3
        lg.legend.spacing = 1

    grid = gridplot([[p, lg]])

    p.frame_width = int(p.plot_width * 0.8)

    # save(grid)

    # show(grid)
    if save_svg:
        path = Path(save_svg)
        cache_dir = (path.parent / '._cache')
        cache_dir.mkdir(exist_ok=True)
        heat_svg = cache_dir / 'heat_main_tmp.svg'
        legend_svg = cache_dir / 'heat_legend_tmp.svg'

        p.output_backend = 'svg'
        export_svgs(p, filename=heat_svg)
        lg.output_backend = 'svg'
        export_svgs(lg, filename=legend_svg)

        Figure(f"{pw + 200}px", f"{ph}px",
               SVG(heat_svg),
               SVG(legend_svg).move(pw, 0),
               ).save(save_svg)
        # .tile(2, 1)

        shutil.rmtree(cache_dir)
    # show(row(p,lg))


