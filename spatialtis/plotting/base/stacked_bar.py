from typing import Any, Dict, List, Optional, Union

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from bokeh.models import ColumnDataSource, FactorRange, Legend
from bokeh.plotting import figure

from spatialtis.plotting.abc import BokehMixin, MatplotlibMixin

from ...utils import doc
from .palette import get_colors


@doc
class stacked_bar_interactive(BokehMixin):
    """Stacked bar plot, Bokeh

    Args:
        data: {data_df}
        groupby: {groupby}
        stacked_types: The level used to stacked
        percentage: To scale the data to percentage
        sort_by: Sort by which value
        ascending: Ascending sort
        group_order: {group_order}
        direction: {direction}
        **plot_options: {plot_options}

    """

    def __init__(
        self,
        data: pd.DataFrame,
        groupby: List,
        stacked_types: List[str],
        percentage: bool = True,
        sort_by: Union[List[str], str, None] = None,
        ascending: bool = True,
        group_order: Optional[Dict[str, List]] = None,
        direction: Union[str] = "vertical",
        **plot_options,
    ):
        super().__init__(**plot_options)
        self.data = data
        self.groupby = groupby
        self.group_order = group_order

        gl = len(groupby)
        if gl > 3:
            raise ValueError("Only support 3 levels depth categorical data for bokeh")

        if group_order is not None:
            data = self.reorder(data)
        data = data.groupby(groupby, sort=False).sum()
        plot_data = data[stacked_types].reset_index(drop=True)

        if percentage:
            plot_data = plot_data.div(plot_data.sum(axis=1), axis=0) * 100
            plot_data.fillna(0)
        if sort_by:
            plot_data = plot_data.sort_values(sort_by, ascending=ascending)
        reg = plot_data.to_dict(orient="list")
        factors = list(data.index)
        reg["factors"] = factors

        types_count = len(stacked_types)

        # config for figure
        figure_config: Dict[str, Any] = dict(
            tools="save,hover",
            toolbar_location=None,
            tooltips="$name @$name%" if percentage else "$name @$name",
        )

        if self.size is None:
            if direction == "vertical":
                figure_config["plot_height"] = 400
            elif direction == "horizontal":
                figure_config["plot_width"] = 400
        else:
            figure_config["plot_height"] = self.size[0]
            figure_config["plot_width"] = self.size[1]

        # set colors
        default_palette = ["Category20", "Spectral"]
        if self.palette is None:
            palette = default_palette
        colors = get_colors(types_count, palette)

        franger = FactorRange(
            *factors, group_padding=0, factor_padding=-0.45, subgroup_padding=-0.35
        )
        if direction == "vertical":
            figure_config["x_range"] = factors if gl == 1 else franger
        else:
            figure_config["y_range"] = factors if gl == 1 else franger

        # do the plotting
        source = ColumnDataSource(data=reg)
        p = figure(**figure_config)

        # some beautify setting

        if direction == "vertical":
            b = p.vbar_stack(
                stacked_types,
                x="factors",
                width=0.5,
                alpha=0.7,
                color=colors,
                source=source,
            )
            p.y_range.start = 0
            p.xaxis.major_label_orientation = 1
            p.xgrid.grid_line_color = None
        else:
            b = p.hbar_stack(
                stacked_types,
                y="factors",
                height=0.5,
                alpha=0.7,
                color=colors,
                source=source,
            )
            p.x_range.start = 0
            p.yaxis.major_label_orientation = 1
            p.ygrid.grid_line_color = None

        # set legend
        legend = Legend(
            items=[(t, [b[i]]) for i, t in enumerate(stacked_types)],
            location="center_right",
        )
        p.add_layout(legend, "right")
        p.hover.point_policy = "follow_mouse"
        p.xaxis.axis_label = self.xaxis_title
        p.yaxis.axis_label = self.yaxis_title

        self.plot = p
        self.set_up()


@doc
class stacked_bar_static(MatplotlibMixin):
    """Stacked bar plot, Matplotlib

    Args:
        data: {data_df}
        groupby: {groupby}
        stacked_types: The level used to stacked
        percentage: To scale the data to percentage
        sort_by: Sort by which value
        ascending: Ascending sort
        group_order: {group_order}
        direction: {direction}
        **plot_options: {plot_options}

    """

    def __init__(
        self,
        data: pd.DataFrame,
        groupby: List,
        stacked_types: List[str],
        percentage: bool = True,
        sort_by: Union[List[str], str, None] = None,
        ascending: bool = True,
        group_order: Optional[Dict[str, List]] = None,
        direction: str = "vertical",
        **plot_options,
    ):
        super().__init__(**plot_options)

        gl = len(groupby)
        if gl > 2:
            raise ValueError(
                "Only support 2 levels depth categorical data for matplotlib"
            )

        data = data[groupby + stacked_types]
        self.data = data
        self.groupby = groupby
        self.group_order = group_order
        plot_data = data.groupby(groupby, sort=False).sum().reset_index()

        if sort_by is not None:
            plot_data = plot_data.sort_values(sort_by, ascending=ascending)
        if group_order is not None:
            plot_data = self.reorder(plot_data)

        c = plot_data[stacked_types].T
        if percentage:
            c = c / c.sum() * 100
        c = c.T
        c_cols = c.columns
        for ix, col in enumerate(c_cols):
            if ix > 0:
                c[col] += c[c_cols[ix - 1]]
        plot_data[stacked_types] = c

        default_palette = ["Category20", "Spectral"]
        if self.palette is None:
            self.palette = default_palette
        colors = get_colors(len(stacked_types), self.palette)

        level1 = groupby[0]
        if gl > 1:
            level2 = groupby[1]
            groups = plot_data.groupby(level1, sort=False)
            self.fig, self.axes = plt.subplots(1, len(groups), figsize=self.size)
            if direction == "vertical":
                for ix, ((name, g), ax) in enumerate(zip(groups, self.axes.flatten())):
                    for c, co in zip(stacked_types[::-1], colors):
                        sg = sns.barplot(
                            x=level2,
                            y=c,
                            data=g,
                            label=c,
                            ax=ax,
                            color=co,
                            order=g[level2],
                        )
                    sg.set_xticklabels(
                        ax.get_xticklabels(),
                        rotation=self.xtickslabel_rotation,
                        ha=self.xtickslabel_loc,
                    )
                    if level2 == level1:
                        title_name = None
                    else:
                        title_name = name
                    sg.set(ylabel="", xlabel="", title=title_name)
                    if ix != 0:
                        sg.set(yticklabels=[])
                        sg.tick_params(left=False)
            else:
                for ix, ((name, g), ax) in enumerate(zip(groups, self.axes.flatten())):
                    for c, co in zip(stacked_types[::-1], colors):
                        sg = sns.barplot(
                            y=level2,
                            x=c,
                            data=g,
                            label=c,
                            ax=ax,
                            color=co,
                            order=g[level2],
                        )
                    sg.set_yticklabels(
                        ax.get_yticklabels(),
                        rotation=self.ytickslabel_rotation,
                        ha=self.ytickslabel_loc,
                    )
                    if level2 == level1:
                        title_name = None
                    else:
                        title_name = name
                    sg.set(ylabel="", xlabel="", title=title_name)
                    if ix != 0:
                        sg.set(xticklabels=[])
                        sg.tick_params(left=False)
        else:
            self.fig, self.ax = plt.subplots(figsize=self.size)
            if direction == "vertical":
                for c, co in zip(stacked_types[::-1], colors):
                    sb = sns.barplot(
                        x=level1,
                        y=c,
                        data=plot_data,
                        label=c,
                        ax=self.ax,
                        color=co,
                        order=plot_data[level1],
                    )
                sb.set_xticklabels(
                    self.ax.get_xticklabels(),
                    rotation=self.xtickslabel_rotation,
                    ha=self.xtickslabel_loc,
                )
                sb.set(ylabel="", xlabel="")
            else:
                for c, co in zip(stacked_types[::-1], colors):
                    sb = sns.barplot(
                        y=level1,
                        x=c,
                        data=plot_data,
                        label=c,
                        ax=self.ax,
                        color=co,
                        order=plot_data[level1],
                    )
                sb.set_yticklabels(
                    self.ax.get_yticklabels(),
                    rotation=self.ytickslabel_rotation,
                    ha=self.xtickslabel_loc,
                )
                sb.set(ylabel="", xlabel="")

        self.fig.text(0.5, -0.05, self.xaxis_title, ha="center")
        self.fig.text(-0.05, 0.5, self.yaxis_title, va="center", rotation="vertical")
        plt.tight_layout()
        plt.legend(
            bbox_to_anchor=(1.1, 0), borderaxespad=0, loc="lower left", frameon=False
        )
        self.set_up()
