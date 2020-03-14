from typing import Optional, Sequence, Union

import numpy as np
import pandas as pd
from bokeh.io import export_svgs, output_file, output_notebook, save
from bokeh.models import FactorRange, Legend, LegendItem
from bokeh.plotting import figure, show
from scipy.stats import gaussian_kde as kde

from spatialtis import CONFIG

from .palette import get_colors


class Status(object):
    gl: Optional[float] = None
    factors: Optional[Sequence] = None
    q1: pd.DataFrame = None
    q2: pd.DataFrame = None
    q3: pd.DataFrame = None
    upper: pd.DataFrame = None
    lower: pd.DataFrame = None
    splits: Optional[list] = None


s = Status()


def _quantilefy(groups, target_col):
    s.q1 = groups.quantile(q=0.25)
    s.q2 = groups.quantile(q=0.5)
    s.q3 = groups.quantile(q=0.75)
    iqr = s.q3 - s.q1
    s.upper = s.q3 + 1.5 * iqr
    s.lower = s.q1 - 1.5 * iqr
    qmin = groups.quantile(q=0.00)
    qmax = groups.quantile(q=1.00)
    s.upper[target_col] = [
        min([x, y])
        for (x, y) in zip(list(qmax.loc[:, target_col]), s.upper[target_col])
    ]
    s.lower[target_col] = [
        max([x, y])
        for (x, y) in zip(list(qmin.loc[:, target_col]), s.lower[target_col])
    ]


def _kde_points(raw_data, plot_points):
    kernel = kde(raw_data)
    return kernel.evaluate(plot_points)


def _violin_patches(
    groups: pd.core.groupby.DataFrameGroupBy, target_col: str, side: str = "both",
):
    violins = list()
    padding = 0.5
    for i, (_, g) in enumerate(groups):
        data = g[target_col].to_numpy()
        if not len(data) > 1:
            raise ValueError("Can't plot with only one point for each group")
        split_interval = round(len(data) / 4)
        split_interval = split_interval if split_interval >= 100 else 100
        points = np.linspace(np.min(data), np.max(data), split_interval)
        curvepoint = _kde_points(data, points)
        # set endpoints to zero to close violin patches
        norm_curve = (
            (curvepoint - np.min(curvepoint))
            / (np.max(curvepoint) - np.min(curvepoint))
        ) * 0.3
        norm_curve[0] = 0
        norm_curve[-1] = 0
        if side == "both":
            violins.append(
                [
                    np.hstack((-norm_curve + i + padding, norm_curve + i + padding,)),
                    np.hstack((points, points,)),
                ]
            )
        elif side == "+":
            violins.append(
                [np.hstack((norm_curve + i + padding,)), np.hstack((points,)),]
            )
        elif side == "-":
            violins.append(
                [np.hstack((-norm_curve + i + padding,)), np.hstack((points,)),]
            )

    return violins


def _violin_split(violins_a, violins_b, figure_config, colors, direction="vertical"):
    p = figure(**figure_config)

    legends = list()
    for i, (v1, v2) in enumerate(zip(violins_a, violins_b)):
        if direction == "vertical":
            x1, y1 = v1[0], v1[1]
            x2, y2 = v2[0], v2[1]
        elif direction == "horizontal":
            x1, y1 = v1[1], v1[0]
            x2, y2 = v2[1], v2[0]

        b1 = p.patch(x1, y1, fill_color=colors[0], line_color="black", alpha=0.8)
        b2 = p.patch(x2, y2, fill_color=colors[1], line_color="black", alpha=0.8)
        if i == 0:
            legends.append(LegendItem(label=s.splits[0], renderers=[b1]))
            legends.append(LegendItem(label=s.splits[1], renderers=[b2]))

    p.add_layout(Legend(items=legends, location="center_right"), "right")

    return p


def _violin_main(violins, target_col, figure_config, colors, direction="vertical"):
    p = figure(**figure_config)

    for i, v in enumerate(violins):
        if direction == "vertical":
            x, y = v[0], v[1]
        elif direction == "horizontal":
            x, y = v[1], v[0]

        if s.gl > 1:
            p.patch(
                x,
                y,
                fill_color=colors[s.factors[i][-1]],
                line_color=colors[s.factors[i][-1]],
            )
        else:
            p.patch(
                x, y, fill_color=colors[s.factors[i]], line_color=colors[s.factors[i]]
            )

    q23 = s.q2[target_col] - s.q3[target_col]
    q12 = s.q1[target_col] - s.q2[target_col]

    if direction == "vertical":
        p.segment(
            s.factors,
            s.upper[target_col],
            s.factors,
            s.lower[target_col],
            line_color="black",
        )

        p.rect(
            s.factors,
            s.q3[target_col] + q23 / 2,
            0.05,
            q23,
            fill_color="#E08E79",
            line_color="black",
        )
        p.rect(
            s.factors,
            s.q2[target_col] + q12 / 2,
            0.05,
            q12,
            fill_color="#3B8686",
            line_color="black",
        )
    elif direction == "horizontal":
        p.segment(
            s.upper[target_col],
            s.factors,
            s.lower[target_col],
            s.factors,
            line_color="black",
        )

        p.rect(
            s.q3[target_col] + q23 / 2,
            s.factors,
            q23,
            0.05,
            fill_color="#E08E79",
            line_color="black",
        )
        p.rect(
            s.q2[target_col] + q12 / 2,
            s.factors,
            q12,
            0.05,
            fill_color="#3B8686",
            line_color="black",
        )

    return p


def _set_figure_config(title=None, size=None, direction="vertical"):
    # config for figure
    figure_config = dict(tools="save", toolbar_location=None, title=title)
    franger = FactorRange(*s.factors, group_padding=0, subgroup_padding=0)
    if direction == "vertical":
        figure_config["x_range"] = s.factors if s.gl == 1 else franger
    elif direction == "horizontal":
        figure_config["y_range"] = s.factors if s.gl == 1 else franger

    if size is None:
        if direction == "vertical":
            figure_config["plot_height"] = 400
        elif direction == "horizontal":
            figure_config["plot_width"] = 400
    else:
        figure_config["plot_height"] = size[0]
        figure_config["plot_width"] = size[1]

    return figure_config


def _set_colors(palette, mapper=False):
    default_palette = ["Set3"]
    if palette is None:
        palette = default_palette
    if not mapper:
        return get_colors(2, palette)

    if s.gl > 1:
        unique_factor = pd.unique([i[-1] for i in s.factors])
        colors = dict(zip(unique_factor, get_colors(len(unique_factor), palette)))
    else:
        colors = dict(zip(s.factors, get_colors(len(s.factors), palette)))
    return colors


def violin_plot(
    df: pd.DataFrame,
    groupby: Sequence,
    target_col: str,
    split: Optional[str] = None,
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
        raise ValueError(
            f"Unrecognized direction '{direction}', available options are `vertical` and `horizontal`."
        )

    s.gl = len(groupby)

    if isinstance(split, str):
        if split in groupby:
            s.splits = np.unique(dict(zip(df.index.names, df.index.levels))[split])
            if len(s.splits) != 2:
                raise ValueError(
                    f"Can't split more than 2 distinct elements in column '{split}'"
                )
            groupby.remove(split)
            s.gl = s.gl - 1

            if s.gl > 3:
                raise ValueError(
                    "Only support 3 levels depth categorical data, maybe too much group_by elements"
                )

            groups = df.loc[:, [target_col]].groupby(level=split)

            subg_a = groups.get_group(s.splits[0]).groupby(level=groupby)
            violins_a = _violin_patches(subg_a, target_col, side="-")

            subg_b = groups.get_group(s.splits[1]).groupby(level=groupby)
            violins_b = _violin_patches(subg_b, target_col, side="+")

            s.factors = [n for n, _ in subg_a]
            figure_config = _set_figure_config(
                title=title, size=size, direction=direction
            )
            colors = _set_colors(palette)

            p = _violin_split(
                violins_a, violins_b, figure_config, colors, direction=direction
            )

        else:
            raise Exception(f"Index name not exist, '{split}'")
    else:
        if s.gl > 3:
            raise Exception("Only support 3 levels depth categorical data")

        groups = df.loc[:, [target_col]].groupby(level=groupby)
        _quantilefy(groups, target_col)

        violins = _violin_patches(groups, target_col)
        s.factors = [n for n, _ in groups]
        figure_config = _set_figure_config(title=title, size=size, direction=direction)
        colors = _set_colors(palette, mapper=True)

        p = _violin_main(
            violins, target_col, figure_config, colors, direction=direction
        )

    if direction == "vertical":
        p.xgrid.grid_line_color = None
        p.ygrid.grid_line_alpha = 0.7
    elif direction == "horizontal":
        p.ygrid.grid_line_color = None
        p.xgrid.grid_line_alpha = 0.7

    # save something
    if save_html:
        output_file(save_html)
        save(p)

    if save_svg:
        p.output_backend = "svg"
        export_svgs(p, filename=save_svg)

    if (CONFIG.WORKING_ENV is not None) & display:
        output_notebook(hide_banner=True, notebook_type=CONFIG.WORKING_ENV)
        show(p)

    if return_plot:
        return p
