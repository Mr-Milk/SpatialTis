from typing import Optional, Sequence, Union

import numpy as np
import pandas as pd
from bokeh.io import export_svgs, output_file, output_notebook, save, show
from bokeh.models import Legend, LegendItem
from bokeh.plotting import figure

from ..config import WORKING_ENV
from .palette import get_colors


def cell_map(
    df: pd.DataFrame,
    type_col: str,
    selected_types: Optional[Sequence] = None,
    shape_key: Optional[str] = "cell_shape",
    display: bool = True,
    title: Optional[str] = None,
    save_svg: Optional[str] = None,
    save_html: Optional[str] = None,
    palette: Union[Sequence[str], str, None] = None,
    return_plot: bool = False,
):

    groups = df.groupby(type_col)

    default_palette = "Spectral", "Category20"
    if palette is None:
        palette = default_palette
    colors = get_colors(len(groups), *palette)

    tools = "pan,wheel_zoom,box_zoom,reset,hover,save"

    p = figure(
        title=title,
        tools=tools,
        x_axis_location=None,
        y_axis_location=None,
        toolbar_location="above",
        tooltips="@name",
    )

    legends = list()

    def add_patches(name, fill_color=None, fill_alpha=None):
        x = [[c[0] for c in eval(cell)] for cell in data[shape_key]]
        y = [[c[1] for c in eval(cell)] for cell in data[shape_key]]
        plot_data = dict(x=x, y=y, name=[name for i in range(0, len(x))])
        b = p.patches(
            "x",
            "y",
            source=plot_data,
            fill_color=fill_color,
            fill_alpha=fill_alpha,
            line_color="white",
            line_width=0.5,
        )
        legends.append(LegendItem(label=name, renderers=[b]))

    if selected_types is None:
        for ix, (n, data) in enumerate(groups):
            add_patches(n, fill_color=colors[ix], fill_alpha=0.8)
    else:
        for ix, (n, data) in enumerate(groups):
            if n in selected_types:
                add_patches(n, fill_color=colors[ix], fill_alpha=0.8)
            else:
                add_patches("other", fill_color="grey", fill_alpha=0.5)

    p.add_layout(Legend(items=legends, location="center_right"), "right")

    p.grid.grid_line_color = None
    p.hover.point_policy = "follow_mouse"
    p.legend.label_text_font_size = "8pt"
    p.legend.glyph_width = 10
    p.legend.glyph_height = 10
    p.legend.click_policy = "hide"
    p.legend.label_text_baseline = "bottom"
    p.legend.spacing = 1

    # save something
    if save_html:
        output_file(save_html)
        save(p)

    if save_svg:
        p.output_backend = "svg"
        export_svgs(p, filename=save_svg)

    # solve env here
    if (WORKING_ENV is not None) & display:
        output_notebook(hide_banner=True, notebook_type=WORKING_ENV)
        show(p)

    # it will return a bokeh plot instance, allow user to do some modification
    if return_plot:
        return p
