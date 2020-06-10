from typing import Any, Dict, Mapping, Optional, Sequence, Union

import pandas as pd
from bokeh.io import output_notebook
from bokeh.models import ColumnDataSource, FactorRange, Legend
from bokeh.plotting import figure, show

from spatialtis import CONFIG

from ._save import save_bokeh
from .palette import get_colors


def stacked_bar(
    df: pd.DataFrame,
    groupby: Sequence[str],
    percentage: bool = True,
    sort_type: Optional[str] = None,
    ascending: bool = True,
    group_order: Optional[tuple] = None,
    direction: Union[str] = "vertical",
    size: Optional[Sequence[int]] = None,
    title: Optional[str] = None,
    palette: Union[Sequence[str], str, None] = None,
    display: bool = True,
    save: Optional[str] = None,
    return_plot: bool = False,
):
    if direction not in ["vertical", "horizontal"]:
        raise ValueError(f"Unrecognized direction '{direction}'")

    gl = len(groupby)

    if gl > 3:
        print("Only support 3 levels depth categorical data")
        return None

    df = df.groupby(level=groupby).sum()

    if percentage:
        df = df.div(df.sum(axis=1), axis=0) * 100

    if sort_type is not None:
        df = df.sort_values(sort_type, ascending=ascending)

    if group_order is not None:
        mapper = dict()
        reorder_index = []
        na_value = max([len(i) for i in group_order.values()])
        for ix in df.index:
            nums = []
            if isinstance(ix, str):
                nums.append(group_order[groupby[0]].index(ix))
            else:
                for i, e in enumerate(ix):
                    try:
                        nums.append(group_order[groupby[i]].index(e))
                    except:
                        nums.append(na_value + 1)
            try:
                mapper[tuple(nums)].append(ix)
            except:
                mapper[tuple(nums)] = []
                mapper[tuple(nums)].append(ix)

        for i in sorted(list(mapper.keys())):
            for t in mapper[i]:
                reorder_index.append(t)
        df = df.loc[reorder_index, :]

    factors = list()
    for i in list(df.index):
        if isinstance(i, Sequence):
            if not isinstance(i, str):
                factors.append(tuple([str(t) for t in i]))
            else:
                factors.append(i)
        else:
            factors.append(str(i))

    types = list(df.columns)
    reg = df.to_dict(orient="list")

    reg["factors"] = factors

    bar_count = len(reg[list(reg.keys())[0]])
    types_count = len(types)

    # config for figure
    figure_config: Dict[str, Any] = dict(
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
    if save is not None:
        save_bokeh(p, save)

    # solve env here
    if (CONFIG.WORKING_ENV is not None) & display:
        output_notebook(hide_banner=True, notebook_type=CONFIG.WORKING_ENV)
        show(p)

    # it will return a bokeh plot instance, allow user to do some modification
    if return_plot:
        return p
