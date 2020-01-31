from collections import OrderedDict
from typing import Sequence, Union, Optional

import numpy as np
import pandas as pd
from bokeh.plotting import figure, show
from bokeh.models import ColumnDataSource, FactorRange, Legend
from bokeh.io import output_notebook, export_svgs, save

from spatialtis.plotting.palette import get_colors, colorcycle
from spatialtis.config import WORKING_ENV


def regroup_df(df: pd.DataFrame, group_by: Sequence[str]):
    gl = len(group_by)

    col = np.array([np.array(i, dtype=str) for i in df.columns.to_numpy()])

    colname = df.columns.names

    depth = df.columns.nlevels

    types = list(df.index)
    sel = [colname.index(g) for g in group_by]
    unsel = set(range(0, len(colname))) - set(sel)
    selected_col = col

    if gl < depth:
        redf = OrderedDict()
        idx = pd.IndexSlice
        for s in selected_col:
            print(s)
            q = [f"'{i}'" for i in s]
            for u in unsel:
                q.insert(u, ':')
                if gl == 1:
                    exec(f'''redf[{str(s[0])}] = np.sum(df.loc[:,idx[{",".join(q)}]].to_numpy(), axis=1)''')
                else:
                    exec(f'''redf[{tuple(s)}] = np.sum(df.loc[:,idx[{",".join(q)}]].to_numpy(), axis=1)''')
        if gl == 1:
            factors = [str(k) for k in redf.keys()]
        else:
            factors = list(redf.keys())
        rowdata = dict(zip(types, np.array(list(redf.values())).T))
        rowdata['type_key'] = factors

        return rowdata, factors, types
    else:
        factors = [tuple(s) for s in selected_col]
        rowdata = dict(zip(df.index, df.to_numpy()))
        rowdata['type_key'] = factors

        return rowdata, factors, types


def stacked_bar(
        df: pd.DataFrame,
        group_by: Sequence[str],
        percentage: bool = True,
        size: Union[Sequence[int], str] = None,
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

    reg, factors, types = regroup_df(df, group_by)

    if percentage:
        percentage_reg = {k: v / np.sum(v) * 100 for k, v in reg.items()}
        reg = percentage_reg

    bar_count = len(reg[list(reg.keys())[0]])
    types_count = len(types)

    # config for figure
    tools = 'save,hover'
    toolbar_location = 'above'
    tooltips = '$name @$name%' if percentage else '$name @$name'  # add % when use percentage
    if size is None:
        plot_height = 400
        plot_width = 400 + (bar_count - 8) * 50 if (bar_count > 8) else 400  # auto-adjust figure
    else:
        plot_height = size[0]
        plot_width = size[1]

    # config plot
    if color is None & palette is None:
        colors = get_colors(colorcycle('Spectral', 'Category20'), types_count)
    elif color is not None:
        colors = color
    else:
        colors = get_colors(colorcycle(palette), types_count)

    # do the plotting
    source = ColumnDataSource(data=reg)
    if gl == 1:
        p = figure(x_range=factors, plot_height=plot_height, plot_width=plot_width,
                   toolbar_location=toolbar_location, tools=tools, tooltip=tooltips)
    else:
        p = figure(x_range=FactorRange(*factors, group_padding=0, factor_padding=-0.45, subgroup_padding=-0.35),
                   plot_height=plot_height, plot_width=plot_width,
                   toolbar_location=toolbar_location, tools=tools, tooltips=tooltips)

    # some beautify setting
    if bar_direction == 'vertical':
        b = p.vbar_stack(types, x="type_key", width=0.5, alpha=0.5, color=colors, source=source)
        p.y_range.start = 0
        p.xaxis.major_label_orientation = 1
        p.xgrid.grid_line_color = None
    elif bar_direction == 'horizontal':
        b = p.hbar_stack(types, y="type_key", height=0.5, alpha=0.5, color=colors, source=source)
        p.x_range.start = 0
        p.yaxis.major_label_orientation = 1
        p.ygrid.grid_line_color = None

    # set legend
    legend = Legend(items=[(t, [b[i]]) for i, t in enumerate(types)], location='center_right')
    p.add_layout(legend, 'right')

    # save something
    if save_html:
        save(p, save_html)

    if save_svg:
        export_svgs(p, filename=save_svg)

    # solve env here
    if WORKING_ENV is None & notebook is not None:
        output_notebook(hide_banner=True, notebook_type=notebook)
        show(p)
