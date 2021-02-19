from typing import Dict, List, Optional, Sequence

import matplotlib.pyplot as plt
import seaborn as sns
from bokeh.models import Legend, LegendItem
from bokeh.plotting import figure
from matplotlib.collections import PatchCollection
from matplotlib.lines import Line2D
from matplotlib.patches import Polygon
from shapely.geometry import MultiPoint

from spatialtis.plotting.abc import BokehMixin, MatplotlibMixin
from spatialtis.plotting.base.palette import get_colors
from spatialtis.utils import doc


@doc
class cell_map_interactive(BokehMixin):
    """Points or Polygons on 2D, Bokeh

    Args:
        points: Array of 2D points
        shapes: Array of 2D polygons
        selected_types: Only selected types will be highlighted, other will be in grey
        cell_size: The size of the points
        **plot_options: {plot_options}

    """

    def __init__(
        self,
        points: Optional[Dict[str, List]] = None,
        shapes: Optional[Dict[str, List]] = None,
        selected_types: Optional[Sequence] = None,
        cell_size: int = 5,
        **plot_options,
    ):
        super().__init__(**plot_options)
        self.cell_size = cell_size
        if points is not None:
            self.data = points
            add_shape = self.add_circle
        else:
            self.data = shapes
            add_shape = self.add_patches

        cell_types = self.data.keys()

        default_palette = ["Category20", "Spectral"]
        if self.palette is None:
            self.palette = default_palette
        colors = get_colors(len(cell_types), self.palette)

        tools = "pan,wheel_zoom,box_zoom,reset,hover,save"

        figure_config = dict(
            title=self.title,
            tools=tools,
            x_axis_location=None,
            y_axis_location=None,
            toolbar_location="above",
            tooltips="@name",
        )

        if self.size is None:
            figure_config["plot_height"] = 700
        else:
            figure_config["plot_height"] = self.size[0]
            figure_config["plot_width"] = self.size[1]

        self.plot = figure(**figure_config)

        self.legends = list()
        self.legends_name: List[str] = list()

        if selected_types is None:
            selected_types = cell_types
        for color, (n, d) in zip(colors, self.data.items()):
            if n in selected_types:
                add_shape(d, n, fill_color=color, fill_alpha=0.8)
            else:
                add_shape(d, "other", fill_color="grey", fill_alpha=0.5)

        if len(self.legends) >= 14:
            cut = int(len(self.legends) // 2)
            legends1 = self.legends[0:cut]
            legends2 = self.legends[cut::]
            self.plot.add_layout(
                Legend(items=legends1, location="center_right"), "right"
            )
            self.plot.add_layout(
                Legend(items=legends2, location="center_right"), "right"
            )

        else:
            self.plot.add_layout(
                Legend(items=self.legends, location="center_right"), "right"
            )

        self.plot.grid.grid_line_color = None
        self.plot.hover.point_policy = "follow_mouse"
        self.plot.legend.label_text_font_size = "8pt"
        self.plot.legend.glyph_width = 10
        self.plot.legend.glyph_height = 10
        self.plot.legend.click_policy = "hide"
        self.plot.legend.label_text_baseline = "bottom"
        self.plot.legend.spacing = 1
        self.plot.match_aspect = True

        self.set_up()

    def add_patches(self, data, name, fill_color=None, fill_alpha=None):
        x = [[c[0] for c in cell] for cell in data]
        y = [[c[1] for c in cell] for cell in data]
        plot_data = dict(x=x, y=y, name=[name for _ in range(len(x))])
        b = self.plot.patches(
            "x",
            "y",
            source=plot_data,
            fill_color=fill_color,
            fill_alpha=fill_alpha,
            line_color="white",
            line_width=0.5,
        )
        if name not in self.legends_name:
            self.legends_name.append(name)
            self.legends.append(LegendItem(label=name, renderers=[b]))

    def add_circle(self, data, name, fill_color=None, fill_alpha=None):
        x = [c[0] for c in data]
        y = [c[1] for c in data]
        plot_data = dict(x=x, y=y, name=[name for _ in range(len(x))])

        b = self.plot.circle(
            "x",
            "y",
            source=plot_data,
            fill_color=fill_color,
            fill_alpha=fill_alpha,
            line_color="white",
            line_width=0.5,
            size=self.cell_size,
        )
        if name not in self.legends_name:
            self.legends_name.append(name)
            self.legends.append(LegendItem(label=name, renderers=[b]))


@doc
class cell_map_static(MatplotlibMixin):
    """Points or Polygons on 2D, Bokeh

        Args:
            points: Array of 2D points
            shapes: Array of 2D polygons
            selected_types: Only selected types will be highlighted, other will be in grey
            cell_size: The size of the points
            **plot_options: {plot_options}

        """

    def __init__(
        self,
        points: Optional[Dict[str, List]] = None,
        shapes: Optional[Dict[str, List]] = None,
        selected_types: Optional[Sequence] = None,
        cell_size: int = 5,
        **plot_options,
    ):
        super().__init__(**plot_options)
        self.cell_size = cell_size

        if points is not None:
            self.data = points
        else:
            self.data = shapes

        cell_types = self.data.keys()
        if selected_types is None:
            selected_types = cell_types

        default_palette = ["Spectral", "Category20"]
        if self.palette is None:
            self.palette = default_palette

        self.fig, self.ax = plt.subplots()
        colors = get_colors(len(self.data.keys()), self.palette)
        legends = []
        labels = []
        if points is not None:
            hue = []
            x = []
            y = []
            points_colors = []
            for (name, points), color in zip(self.data.items(), colors):
                if name not in selected_types:
                    name = "other"
                    color = "#d3d3d3"
                hue += [name for _ in range(len(points))]
                labels.append(name)
                points_colors.append(color)
                legends.append(
                    Line2D(
                        [0],
                        [0],
                        linestyle="none",
                        marker="o",
                        markersize=5,
                        markerfacecolor=color,
                    )
                )
                x += [p[0] for p in points]
                y += [p[1] for p in points]

            ss = sns.scatterplot(
                x=x, y=y, hue=hue, s=self.cell_size, ax=self.ax, fc=points_colors
            )
            plt.tight_layout()
        else:
            patches = []
            points = []
            patch_colors = []
            for (name, polys), color in zip(self.data.items(), colors):
                if name not in selected_types:
                    name = "other"
                    color = "#d3d3d3"
                for poly in polys:
                    patches.append(Polygon(poly, label=name))
                    patch_colors.append(color)
                    points += poly
                legends.append(
                    Line2D(
                        [0],
                        [0],
                        linestyle="none",
                        marker="o",
                        markersize=5,
                        markerfacecolor=color,
                    )
                )
                labels.append(name)
            p = PatchCollection(patches, alpha=0.8, facecolors=patch_colors)
            xmin, ymin, xmax, ymax = MultiPoint(points).bounds
            plt.xlim(xmin, xmax)
            plt.ylim(ymin, ymax)
            self.ax.add_collection(p)

        self.ax.set_aspect("equal")
        plt.axis("off")
        self.ax.legend(
            legends,
            labels,
            title=self.legend_title,
            ncol=2,
            bbox_to_anchor=(1, 0.05),
            borderaxespad=0,
            loc="lower left",
            frameon=False,
        )
        self.set_up()
