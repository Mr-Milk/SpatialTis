from typing import Optional, Sequence, Union

import numpy as np
import pandas as pd
from bokeh.io import output_notebook
from bokeh.models import FactorRange, Legend, LegendItem
from bokeh.plotting import figure, show
from scipy.stats import gaussian_kde as kde

from spatialtis import CONFIG
from spatialtis.utils import reuse_docstring

from .palette import get_colors
from .save import save_bokeh


def _kde_points(data, points):
    kernel = kde(data)
    return kernel.evaluate(points)


class ViolinMain:
    def __init__(
        self,
        df,
        groupby,
        target,
        direction,
        title,
        xaxis_title,
        yaxis_title,
        size,
        palette,
    ):

        self.gl = len(groupby)
        self.q1 = None
        self.q2 = None
        self.q3 = None
        self.upper = None
        self.lower = None
        self.groups = df.loc[:, [target]].groupby(level=groupby, sort=False)
        self.factors = [n for n, _ in self.groups]
        self.target = target
        self.direction = direction
        self.title = title
        self.xaxis_title = xaxis_title
        self.yaxis_title = yaxis_title
        self.size = size
        self.palette = palette

        self._quantilefy()
        self.figure_config = self._set_figure_config()
        self.colors = self._set_colors()
        self.p = self._wrap()

    def _quantilefy(self):
        self.q1 = self.groups.quantile(q=0.25)
        self.q2 = self.groups.quantile(q=0.5)
        self.q3 = self.groups.quantile(q=0.75)
        iqr = self.q3 - self.q1
        self.upper = self.q3 + 1.5 * iqr
        self.lower = self.q1 - 1.5 * iqr
        qmin = self.groups.quantile(q=0.00)
        qmax = self.groups.quantile(q=1.00)
        self.upper[self.target] = [
            min([x, y])
            for (x, y) in zip(list(qmax.loc[:, self.target]), self.upper[self.target])
        ]
        self.lower[self.target] = [
            max([x, y])
            for (x, y) in zip(list(qmin.loc[:, self.target]), self.lower[self.target])
        ]

    def _violin_patches(self):
        violins = list()
        padding = 0.5
        for i, (_, g) in enumerate(self.groups):
            # subsampling when data point exceed 100000
            if len(g) > 100000:
                g = g.sample(n=100000, axis=0, random_state=0)
            data = g[self.target].to_numpy()
            if not len(data) > 1:
                violins.append(
                    [
                        np.asarray([i + padding, i + padding]),
                        np.asarray([data[0], data[0]]),
                    ]
                )
            else:
                # reduce actual drawing resolution of curve
                subsample_size = round(len(data) / 4)
                if subsample_size < 100:
                    subsample_size = 100
                elif subsample_size > 500:
                    subsample_size = 500
                points = np.linspace(np.min(data), np.max(data), subsample_size)
                curvepoint = _kde_points(data, points)
                # set endpoints to zero to close violin patches
                norm_curve = (
                    (curvepoint - np.min(curvepoint))
                    / (np.max(curvepoint) - np.min(curvepoint))
                ) * 0.3
                norm_curve[0] = 0
                norm_curve[-1] = 0
                violins.append(
                    [
                        np.hstack(
                            (-norm_curve + i + padding, norm_curve + i + padding,)
                        ),
                        np.hstack((points, points,)),
                    ]
                )

        return violins

    def _wrap(self):
        p = figure(**self.figure_config)

        violins = self._violin_patches()

        for i, v in enumerate(violins):
            if self.direction == "vertical":
                x, y = v[0], v[1]
            elif self.direction == "horizontal":
                x, y = v[1], v[0]

            if self.gl > 1:
                p.patch(
                    x,
                    y,
                    fill_color=self.colors[self.factors[i][-1]],
                    line_color=self.colors[self.factors[i][-1]],
                )
            else:
                p.patch(
                    x,
                    y,
                    fill_color=self.colors[self.factors[0]],
                    line_color=self.colors[self.factors[0]],
                )

        q23 = self.q2[self.target] - self.q3[self.target]
        q12 = self.q1[self.target] - self.q2[self.target]

        if self.direction == "vertical":
            p.segment(
                self.factors,
                self.upper[self.target],
                self.factors,
                self.lower[self.target],
                line_color="black",
            )

            p.rect(
                self.factors,
                self.q3[self.target] + q23 / 2,
                0.05,
                q23,
                fill_color="#E08E79",
                line_color="black",
            )
            p.rect(
                self.factors,
                self.q2[self.target] + q12 / 2,
                0.05,
                q12,
                fill_color="#3B8686",
                line_color="black",
            )
        elif self.direction == "horizontal":
            p.segment(
                self.upper[self.target],
                self.factors,
                self.lower[self.target],
                self.factors,
                line_color="black",
            )

            p.rect(
                self.q3[self.target] + q23 / 2,
                self.factors,
                q23,
                0.05,
                fill_color="#E08E79",
                line_color="black",
            )
            p.rect(
                self.q2[self.target] + q12 / 2,
                self.factors,
                q12,
                0.05,
                fill_color="#3B8686",
                line_color="black",
            )

        if self.direction == "vertical":
            p.xgrid.grid_line_color = None
            p.xaxis.major_label_orientation = 1
            p.ygrid.grid_line_alpha = 0.7
        elif self.direction == "horizontal":
            p.ygrid.grid_line_color = None
            p.yaxis.major_label_orientation = 1
            p.xgrid.grid_line_alpha = 0.7

        p.xaxis.axis_label = self.xaxis_title
        p.yaxis.axis_label = self.yaxis_title

        return p

    def _set_figure_config(self):
        # config for figure
        figure_config = dict(tools="save", toolbar_location=None, title=self.title)
        franger = FactorRange(*self.factors, group_padding=0, subgroup_padding=0)
        if self.direction == "vertical":
            figure_config["x_range"] = self.factors if self.gl == 1 else franger
        elif self.direction == "horizontal":
            figure_config["y_range"] = self.factors if self.gl == 1 else franger

        if self.size is None:
            if self.direction == "vertical":
                figure_config["plot_height"] = 400
            elif self.direction == "horizontal":
                figure_config["plot_width"] = 400
        else:
            figure_config["plot_height"] = self.size[0]
            figure_config["plot_width"] = self.size[1]

        return figure_config

    def _set_colors(self):
        default_palette = ["Set3"]
        if self.palette is None:
            self.palette = default_palette

        if self.gl > 1:
            unique_factor = pd.unique([i[-1] for i in self.factors])
            colors = dict(
                zip(unique_factor, get_colors(len(unique_factor), self.palette))
            )
        else:
            colors = dict(
                zip(self.factors, get_colors(len(self.factors), self.palette))
            )
        return colors


@reuse_docstring()
def violin_plot(
    df: pd.DataFrame,
    groupby: Sequence,
    target_key: str,
    group_order: Optional[dict] = None,
    direction: Union[str] = "vertical",
    display: Optional[bool] = None,
    title: Optional[str] = None,
    xaxis_title: Optional[str] = None,
    yaxis_title: Optional[str] = None,
    size: Optional[Sequence[int]] = None,
    save: Optional[str] = None,
    palette: Union[Sequence[str], str, None] = None,
    return_plot: bool = False,
):
    """(bokeh) Violin plot

    Args:
        df: Input data
        groupby: {groupby}
        target_key: The key stores the data for plotting in `AnnData.obs`
        group_order: {group_order}
        direction: {direction}
        palette: {palette}
        xaxis_title: {xaxis_title}
        yaxis_title: {yaxis_title}
        size: {size}
        display: {display}
        title: {title}
        save: {save}
        return_plot:  {return_plot}

    """
    if direction not in ["vertical", "horizontal"]:
        raise ValueError(
            f"Unrecognized direction '{direction}', available options are `vertical` and `horizontal`."
        )

    if len(groupby) > 3:
        raise Exception("Only support 3 levels depth categorical data")

    if group_order is not None:
        for level, order in group_order.items():
            df = df.reindex(index=order, level=level)

    plot = ViolinMain(
        df,
        groupby,
        target_key,
        direction,
        title,
        xaxis_title,
        yaxis_title,
        size,
        palette,
    )

    # save something
    if save is not None:
        save_bokeh(plot.p, save)

    if display is None:
        if CONFIG.WORKING_ENV is None:
            display = False
        else:
            display = True
    if display:
        show(plot.p)

    if return_plot:
        return plot.p
