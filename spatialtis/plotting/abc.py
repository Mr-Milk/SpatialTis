from dataclasses import dataclass
from itertools import product
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyecharts.options as opts
import seaborn as sns
from bokeh.io import export_png, export_svgs, output_file, save, show
from matplotlib.axes import Axes, SubplotBase
from matplotlib.figure import Figure
from pyecharts.render import make_snapshot
from snapshot_phantomjs import snapshot

from spatialtis.config import CONFIG
from spatialtis.typing import Number
from spatialtis.utils import doc


@doc
@dataclass
class PlotBase:
    """Base class for all plotting instance

    Noticed that following arguments/attributes might not works for all the plotting instances

    Args:
        plot: The plot instance, only works for Bokeh and Pyecharts plot.
        data: The data used for plotting
        groupby: {groupby}
        title: Major title of the plot
        xaxis_title: The title of X axis
        yaxis_title: The title of Y axis
        legend_title: The title of legend
        xtickslabel_rotation: Degree to rotate X-axis's label
        ytickslabel_rotation: Degree to rotate Y-axis's label
        palette: Either a series of color in hex, or
            `name of palettes <https://docs.bokeh.org/en/latest/docs/reference/palettes.html>`_,
        saved_name: Default file name to save the plot
        save_path: Directory to save the plot
        size: Size of the plot
        display: Whether to show the plot


    """

    # plot itself
    plot: Optional[Any] = None

    # data
    data: Union[np.array, pd.DataFrame, List, Dict, None] = None
    groupby: Union[str, List[str], None] = None
    # group_order: Optional[Dict[str, List[str]]] = None
    # direction: Optional[str] = None

    # text control
    title: Optional[str] = None
    xaxis_title: Optional[str] = None
    yaxis_title: Optional[str] = None
    legend_title: Optional[str] = None

    # ticks control
    # ticks: Optional[bool] = None
    xtickslabel_rotation: int = 0
    ytickslabel_rotation: int = 0

    # color control
    palette: Union[str, List[str], None] = None

    # save
    saved_name: Optional[str] = None
    save_path: Union[str, Path, None] = None
    size: Optional[Tuple[Number, Number]] = None

    # render
    display: bool = True

    def __repr__(self):
        return ""

    def preset(self):
        if CONFIG.WORKING_ENV is None:
            self.display = False
        else:
            self.display = True

    def save(self):
        pass

    def set_up(self):
        pass

    def reorder(self, data):
        levels = []
        orders = []
        for level, order in self.group_order.items():
            levels.append(level)
            orders.append(order)

        ix = [i for i in product(*orders)] if len(orders) > 1 else orders[0]
        orders = dict(zip(ix, range(len(ix))))
        return data.sort_values(levels, key=lambda x: x.map(orders)).reset_index(
            drop=True
        )


class MatplotlibMixin(PlotBase):
    """Mixin for Matplotlib

    Attributes:
        fig: Matplotlib.Figure
        ax: Matplotlib.Subplots or Axes
        axes: Array of ax

    """

    fig: Union[Figure, sns.matrix.ClusterGrid]
    ax: Union[SubplotBase, Axes]
    axes: np.ndarray

    def __init__(self, **plot_options):
        self.preset()
        self.xtickslabel_loc = "center"
        self.ytickslabel_loc = "center"
        super().__init__(**plot_options)

    def set_up(self):
        """Handle display and save"""
        if self.title is not None:
            self.fig.suptitle(self.title)
        if self.display:
            plt.show()
        else:
            plt.close()

        if CONFIG.AUTO_SAVE:
            self.save_path = CONFIG.SAVE_PATH / f"{self.saved_name}.png"
        if self.save_path is not None:
            self.save()

    def save(self, path=None, **save_options):
        """To save a plot

        Args:
            path: The directory and file name to save the plot
            **save_options: Pass to matplotlib's
                `savefig <https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.savefig.html>`_

        """
        if path is not None:
            save_path = Path(path)
        elif self.save_path is not None:
            save_path = Path(self.save_path)
        else:
            save_path = CONFIG.SAVE_PATH
        inbuilt_save_options = dict(dpi=300, bbox_inches="tight")
        for k, v in save_options.items():
            inbuilt_save_options[k] = v
        self.fig.savefig(save_path, **inbuilt_save_options)

    def render(self):
        plt.show()


class BokehMixin(PlotBase):
    """Mixin for Bokeh"""

    def __init__(self, **plot_options):
        self.preset()
        super().__init__(**plot_options)

    def render(self):
        show(self.plot)

    def set_up(self):
        """Handle title, display and save"""
        if self.title is not None:
            self.plot.title.text = self.title
        if self.display:
            show(self.plot)

        if CONFIG.AUTO_SAVE:
            self.save_path = CONFIG.SAVE_PATH / f"{self.saved_name}.html"
        if self.save_path is not None:
            self.save()

    def save(self, path=None):
        """To save a plot

        Better options is to save a html file and use screen capture
        to get a static image

        Args:
            path: The directory and file name to save the plot


        """
        if path is not None:
            save_path = Path(path)
        elif self.save_path is not None:
            save_path = Path(self.save_path)
        else:
            save_path = CONFIG.SAVE_PATH
        file_ext = save_path.suffix[1:]
        save_path = str(save_path)

        plot = self.plot
        plot.background_fill_color = None
        plot.border_fill_color = None

        if file_ext not in ["html", "svg", "png"]:
            raise NotImplementedError(
                "Current supported formats are: svg, html and png"
            )
        elif file_ext == "html":
            output_file(save_path)
            save(plot)
        elif file_ext == "svg":
            plot.output_backend = "svg"
            export_svgs(plot, filename=save_path)
        else:
            export_png(plot, filename=save_path)


class PyechartsMixin(PlotBase):
    """Mixin for Pyecharts

    Args:
        theme: Please go to `theme <https://pyecharts.org/#/zh-cn/themes>`_ for details, (Default: "white")
        renderer: "canvas" or "svg", (Default: "canvas")

    """

    def __init__(self, **plot_options):
        self.theme = "white"
        self.renderer = "canvas"
        self.preset()
        super().__init__(**plot_options)
        if self.size is None:
            self.size = (800, 500)

    def render(self):
        """To render a pyecharts plot"""
        self.plot.load_javascript()
        return self.plot.render_notebook()

    def set_up(self):
        """Handle title, display and save"""
        if self.title is not None:
            self.plot.set_global_opts(title_opts=opts.TitleOpts(self.title),)
        if self.display:
            self.plot.load_javascript()

        if CONFIG.AUTO_SAVE:
            self.save_path = CONFIG.SAVE_PATH / f"{self.saved_name}.html"
        if self.save_path is not None:
            self.save()

    def save(self, path=None, **save_options):
        """To save a plot

        Better options is to save a html file and use screen capture
        to get a static image

        Args:
            path: The directory and file name to save the plot
            **save_options: Pass to `make_snapshot <https://pyecharts.org/#/zh-cn/render_images?id=api>`_

        """
        if path is not None:
            save_path = Path(path)
        elif self.save_path is not None:
            save_path = Path(self.save_path)
        else:
            save_path = CONFIG.SAVE_PATH
        file_ext = save_path.suffix[1:]
        self.plot.renderer = "canvas"

        inbuilt_save_options = dict(delay=2, pixel_ratio=5)
        for k, v in save_options.items():
            inbuilt_save_options[k] = v

        if file_ext not in ["html", "png"]:
            raise NotImplementedError("Current supported formats are: png and html")
        elif file_ext == "html":
            self.plot.render(save_path)
        else:
            make_snapshot(
                snapshot, self.plot.render(), str(save_path), **save_options,
            )
