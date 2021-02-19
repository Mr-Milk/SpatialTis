from typing import Dict, List, Optional

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

from ...utils import doc
from ..abc import MatplotlibMixin


@doc
class violin_static(MatplotlibMixin):
    """Violin plot, Matplotlib

    Args:
        data: {data_df}
        groupby: {groupby}
        target: The level of value
        hue: The level to split value
        group_order: {group_order}
        direction: {direction}
        **plot_options: {plot_options}

    """

    def __init__(
        self,
        data: pd.DataFrame,
        groupby: List,
        target: str,
        hue: Optional[str] = None,
        group_order: Optional[Dict[str, List]] = None,
        direction: Optional[str] = "vertical",
        **plot_options,
    ):
        super().__init__(**plot_options)
        gl = len(groupby)
        if gl > 2:
            raise ValueError(
                "Only support 2 levels depth categorical data for matplotlib"
            )
        self.data = data
        if group_order is not None:
            self.group_order = group_order
            self.data = self.reorder(self.data)
        self.target = target
        self.direction = direction

        level1 = groupby[0]
        if gl > 1:
            level2 = groupby[1]
            groups = self.data.groupby(level1, sort=False)
            self.fig, self.axes = plt.subplots(1, len(groups), figsize=self.size)
            if direction == "vertical":
                for ix, ((name, g), ax) in enumerate(zip(groups, self.axes.flatten())):
                    sg = sns.violinplot(x=level2, y=self.target, data=g, hue=hue, ax=ax)
                    ax.get_legend().remove()
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
                    sg = sns.violinplot(y=level2, x=self.target, data=g, hue=hue, ax=ax)
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
            plt.legend(
                bbox_to_anchor=(1.1, 0),
                borderaxespad=0,
                loc="lower left",
                frameon=False,
            )
        else:
            self.fig, self.ax = plt.subplots(figsize=self.size)
            if direction == "vertical":
                sb = sns.violinplot(x=level1, y=self.target, data=self.data, ax=self.ax)
                sb.set_xticklabels(
                    self.ax.get_xticklabels(),
                    rotation=self.xtickslabel_rotation,
                    ha=self.xtickslabel_loc,
                )
                sb.set(ylabel="", xlabel="")
            else:
                sb = sns.violinplot(y=level1, x=self.target, data=self.data, ax=self.ax)
                sb.set_yticklabels(
                    self.ax.get_yticklabels(),
                    rotation=self.ytickslabel_rotation,
                    ha=self.ytickslabel_loc,
                )
                sb.set(ylabel="", xlabel="")
        self.fig.text(0.5, -0.05, self.xaxis_title, ha="center")
        self.fig.text(-0.05, 0.5, self.yaxis_title, va="center", rotation="vertical")
        plt.tight_layout()
        self.set_up()
