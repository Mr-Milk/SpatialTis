from typing import Optional, Sequence, Union

import pandas as pd
from bokeh.io import export_svgs, output_file, output_notebook, save
from bokeh.models import ColumnDataSource, FactorRange, Legend
from bokeh.plotting import figure, show

from spatialtis import CONFIG
from spatialtis.plotting.palette import get_colors


def _regroup_df(
    df: pd.DataFrame, group_by: Sequence[str], percentage: bool = True,
):
    groups = df.groupby(level=group_by).sum()

    if percentage:
        groups = groups.div(groups.sum(axis=1), axis=0) * 100

    return groups


def stacked_bar(
    df: pd.DataFrame,
    group_by: Sequence[str],
    percentage: bool = True,
    direction: Union[str] = "vertical",
    display: bool = True,
    title: Optional[str] = None,
    size: Optional[Sequence[int]] = None,
    save_svg: Optional[str] = None,
    save_html: Optional[str] = None,
    palette: Union[Sequence[str], str, None] = None,
    return_plot: bool = False,
):
    if direction not in ["vertical", "horizontal"]:
        raise ValueError(f"Unrecognized direction '{direction}'")

    gl = len(group_by)

    if gl > 3:
        print("Only support 3 levels depth categorical data")
        return None

    df = _regroup_df(df, group_by, percentage=percentage)

    factors = list(df.index)
    types = list(df.columns)
    reg = df.to_dict(orient="list")

    reg["factors"] = factors

    bar_count = len(reg[list(reg.keys())[0]])
    types_count = len(types)

    # config for figure
    figure_config = dict(
        tools="save,hover",
        toolbar_location=None,
        tooltips="$name @$name%" if percentage else "$name @$name",
        title=title,
    )

    if size is None:
        if direction == "vertical":
            figure_config["plot_height"] = 400
        elif direction == "horizontal":
            figure_config["plot_width"] = 400
    else:
        figure_config["plot_height"] = size[0]
        figure_config["plot_width"] = size[1]

    # set colors
    default_palette = ["Spectral", "Category20"]
    if palette is None:
        palette = default_palette
    colors = get_colors(types_count, palette)

    franger = FactorRange(
        *factors, group_padding=0, factor_padding=-0.45, subgroup_padding=-0.35
    )
    if direction == "vertical":
        figure_config["x_range"] = factors if gl == 1 else franger
    elif direction == "horizontal":
        figure_config["y_range"] = factors if gl == 1 else franger

    # do the plotting
    source = ColumnDataSource(data=reg)
    p = figure(**figure_config)

    # some beautify setting

    if direction == "vertical":
        b = p.vbar_stack(
            types, x="factors", width=0.5, alpha=0.7, color=colors, source=source
        )
        p.y_range.start = 0
        p.xaxis.major_label_orientation = 1
        p.xgrid.grid_line_color = None
    elif direction == "horizontal":
        b = p.hbar_stack(
            types, y="factors", height=0.5, alpha=0.7, color=colors, source=source
        )
        p.x_range.start = 0
        p.yaxis.major_label_orientation = 1
        p.ygrid.grid_line_color = None

    # set legend
    legend = Legend(
        items=[(t, [b[i]]) for i, t in enumerate(types)], location="center_right"
    )
    p.add_layout(legend, "right")
    p.hover.point_policy = "follow_mouse"

    # save something
    if save_html:
        output_file(save_html)
        save(p)

    if save_svg:
        p.output_backend = "svg"
        export_svgs(p, filename=save_svg)

    # solve env here
    if (CONFIG.WORKING_ENV is not None) & display:
        output_notebook(hide_banner=True, notebook_type=CONFIG.WORKING_ENV)
        show(p)

    # it will return a bokeh plot instance, allow user to do some modification
    if return_plot:
        return p
