from collections import OrderedDict
from typing import Sequence, Union, Optional

import numpy as np
import pandas as pd
from bokeh.plotting import figure, show
from bokeh.models import ColumnDataSource, FactorRange, Legend
from bokeh.io import output_notebook, export_svgs, save

from spatialtis.plotting.palette import get_colors, colorcycle
from spatialtis.config import WORKING_ENV


def _regroup_df(
        df: pd.DataFrame,
        group_by: Sequence[str],
        percentage: bool = True,
):
    groups = df.groupby(level=group_by).sum()

    if percentage:
        groups = groups.div(groups.sum(axis=1), axis=0) * 100

    return groups


def stacked_bar(
        df: pd.DataFrame,
        group_by: Sequence[str],
        percentage: bool = True,
        size: Union[Sequence[int], str] = None,
        title: Optional[str] = None,
        notebook: Optional[str] = None,
        save_svg: Optional[str] = None,
        save_html: Optional[str] = None,
        bar_direction: Union[str] = 'horizontal',
        color: Optional[str] = None,
        palette: Optional[str] = None,
):
    gl = len(group_by)

    if gl > 3:
        print('Only support 3 levels depth categorical data')
        return None

    df = _regroup_df(df, group_by, percentage=percentage)

    factors = list(df.index)
    types = list(df.columns)
    reg = df.to_dict(orient='list')

    reg['factors'] = factors

    bar_count = len(reg[list(reg.keys())[0]])
    types_count = len(types)

    # config for figure
    figure_config = dict(
        tools='save,hover',
        toolbar_location='above',
        tooltips='$name @$name%' if percentage else '$name @$name',
    )

    if size is None:
        if bar_direction == 'vertical':
            figure_config['plot_height'] = 400
            figure_config['plot_width'] = 400 + (bar_count - 8) * 50 if (bar_count > 8) else 400  # auto-adjust figure
        elif bar_direction == 'horizontal':
            figure_config['plot_width'] = 400
            figure_config['plot_height'] = 400 + (bar_count - 8) * 50 if (bar_count > 8) else 400  # auto-adjust figure
    else:
        figure_config['plot_height'] = size[0]
        figure_config['plot_width'] = size[1]

    if title is not None:
        figure_config['title'] = title

    # config plot
    if (color is None) & (palette is None):
        colors = get_colors(colorcycle('Spectral', 'Category20'), types_count)
    elif color is not None:
        colors = color
    else:
        colors = get_colors(colorcycle(palette), types_count)

    franger = FactorRange(*factors, group_padding=0, factor_padding=-0.45, subgroup_padding=-0.35)
    if bar_direction == 'vertical':
        figure_config['x_range'] = factors if gl == 1 else franger
    elif bar_direction == 'horizontal':
        figure_config['y_range'] = factors if gl == 1 else franger

    # do the plotting
    source = ColumnDataSource(data=reg)
    p = figure(**figure_config)

    # some beautify setting

    if bar_direction == 'vertical':
        b = p.vbar_stack(types, x="factors", width=0.5, alpha=0.7, color=colors, source=source)
        p.y_range.start = 0
        p.xaxis.major_label_orientation = 1
        p.xgrid.grid_line_color = None
    elif bar_direction == 'horizontal':
        b = p.hbar_stack(types, y="factors", height=0.5, alpha=0.7, color=colors, source=source)
        p.x_range.start = 0
        p.yaxis.major_label_orientation = 1
        p.ygrid.grid_line_color = None

    # set legend
    legend = Legend(items=[(t, [b[i]]) for i, t in enumerate(types)], location='center_right')
    p.add_layout(legend, 'right')
    p.hover.point_policy = "follow_mouse"

    # save something
    if save_html:
        save(p, save_html)

    if save_svg:
        export_svgs(p, filename=save_svg)

    # solve env here
    if (WORKING_ENV is None) & (notebook is not None):
        output_notebook(hide_banner=True, notebook_type=notebook)
        show(p)
    elif (WORKING_ENV is not None) & (notebook is None):
        show(p)
