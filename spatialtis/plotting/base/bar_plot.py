import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

from spatialtis.plotting.abc import MatplotlibMixin
from spatialtis.utils import doc


@doc
class bar_static(MatplotlibMixin):
    """Bar plot, Matplotlib

    Args:
        data: {data_df}
        x: The level of value for x axis
        y: The level of value y axis
        sort: To sort the data
        direction: {direction}
        **plot_options: {plot_options}

    """

    def __init__(
        self,
        data: pd.DataFrame,
        x: str,
        y: str,
        sort: bool = True,
        direction: str = "vertical",
        **plot_options,
    ):
        super().__init__(**plot_options)
        if sort:
            data = data.sort_values(x)
        self.fig, self.ax = plt.subplots(figsize=self.size)
        if direction == "vertical":
            sns.barplot(data=data, x=x, y=y, ax=self.ax, ci=None)
        else:
            sns.barplot(data=data, y=x, x=y, ax=self.ax, ci=None)
        plt.tight_layout()
        self.set_up()
