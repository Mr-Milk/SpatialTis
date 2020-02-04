import numpy as np
import pandas as pd

from bokeh.plotting import figure
from bokeh.io import show
import bokeh.palettes as pl
from bokeh.models import Legend, LegendItem
from .palette import colorcycle, get_colors

from typing import Optional, Union, Sequence


def cell_map(df: pd.DataFrame,
             type_col: str,
             selected_types: Optional[Sequence] = None,
             shape_key: Optional[str] = 'cell_shape',
             display: bool = True,
             title: Optional[str] = None,
             size: Optional[Sequence[int]] = None,
             save_svg: Optional[str] = None,
             save_html: Optional[str] = None,
             palette: Union[Sequence[str], str, None] = None,
             ):
    all_types = np.unique(df[type_col])

    colors = get_colors(colorcycle('Spectral', 'Category20'), len(all_types))

    tools = "pan,wheel_zoom,box_zoom,reset,hover,save"

    p = figure(
        title=title, tools=tools,
        x_axis_location=None, y_axis_location=None,
        toolbar_location='above',
        tooltips="@name")

    legends = list()

    def add_patches(name, fill_color=None, fill_alpha=None):
        x = [[c[0] for c in eval(cell)] for cell in data[shape_key]]
        y = [[c[1] for c in eval(cell)] for cell in data[shape_key]]
        plot_data = dict(
            x=x,
            y=y,
            name=[name for i in range(0, len(x))]
        )
        b = p.patches('x', 'y', source=plot_data,
                      fill_color=fill_color,
                      fill_alpha=fill_alpha, line_color="white", line_width=0.5)
        legends.append(LegendItem(label=t, renderers=[b]))

    if selected_types is None:
        for ix, t in enumerate(all_types):
            data = df[df[type_col] == t]
            add_patches(t, fill_color=colors[ix], fill_alpha=0.8)
    else:
        for ix, t in enumerate(selected_types):
            data = df[df[type_col] == t]
            add_patches(t, fill_color=colors[ix], fill_alpha=0.8)
        data = df[~df[type_col].isin(selected_types)]
        add_patches('other', fill_color='grey', fill_alpha=0.5)

    p.add_layout(Legend(items=legends, location='center_right'), 'right')

    p.grid.grid_line_color = None
    p.hover.point_policy = "follow_mouse"
    p.legend.label_text_font_size = '8pt'
    p.legend.glyph_width = 10
    p.legend.glyph_height = 10
    p.legend.click_policy = "hide"
    p.legend.label_text_baseline = 'bottom'
    p.legend.spacing = 1

    show(p)
