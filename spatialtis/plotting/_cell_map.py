from typing import List, Optional, Sequence, Mapping, Union
from pathlib import Path

from anndata import AnnData
from bokeh.io import output_notebook, show
from bokeh.models import Legend, LegendItem
from bokeh.plotting import figure

from spatialtis import CONFIG
from ._save import save_bokeh
from .palette import get_colors


def cell_map(
    adata: AnnData,
    query: Mapping,
    type_col: Optional[str] = None,
    selected_types: Optional[Sequence] = None,
    shape_col: Optional[str] = None,
    size: Optional[Sequence[int]] = None,
    title: Optional[str] = None,
    palette: Union[Sequence[str], str, None] = None,
    display: bool = True,
    save: Union[str, Path, None] = None,
    return_plot: bool = False,
):
    if type_col is None:
        type_col = CONFIG.CELL_TYPE_COL
    if shape_col is None:
        shape_col = CONFIG.SHAPE_COL

    df = adata.obs.query("&".join([f"({k}=='{v}')" for k, v in query.items()]))
    groups = df.groupby(type_col)

    default_palette = ["Spectral", "Category20"]
    if palette is None:
        palette = default_palette
    colors = get_colors(len(groups), palette)

    tools = "pan,wheel_zoom,box_zoom,reset,hover,save"

    figure_config = dict(
        title=title,
        tools=tools,
        x_axis_location=None,
        y_axis_location=None,
        toolbar_location="above",
        tooltips="@name",
    )

    if size is None:
        figure_config["plot_height"] = 700
    else:
        figure_config["plot_height"] = size[0]
        figure_config["plot_width"] = size[1]

    p = figure(**figure_config)

    legends = list()
    legends_name: List[str] = list()

    def add_patches(name, fill_color=None, fill_alpha=None):
        x = [[c[0] for c in eval(cell)] for cell in data[shape_col]]
        y = [[c[1] for c in eval(cell)] for cell in data[shape_col]]
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
        if name not in legends_name:
            legends_name.append(name)
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

    if len(legends) >= 16:
        cut = int(len(legends) // 2)
        legends1 = legends[0:cut]
        legends2 = legends[cut::]
        p.add_layout(Legend(items=legends1, location="center_right"), "right")
        p.add_layout(Legend(items=legends2, location="center_right"), "right")

    else:
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
    if save is not None:
        save_bokeh(p, save)

    # solve env here
    if (CONFIG.WORKING_ENV is not None) & display:
        output_notebook(hide_banner=True, notebook_type=CONFIG.WORKING_ENV)
        show(p)

    # it will return a bokeh plot instance, allow user to do some modification
    if return_plot:
        return p
