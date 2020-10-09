import warnings
from pathlib import Path
from typing import List, Mapping, Optional, Sequence, Union

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
    geom: str = "shape",
    selected_types: Optional[Sequence] = None,
    type_key: Optional[str] = None,
    shape_key: Optional[str] = None,
    centroid_key: Optional[str] = None,
    size: Optional[Sequence[int]] = None,
    title: Optional[str] = None,
    palette: Union[Sequence[str], str, None] = None,
    display: Optional[bool] = None,
    save: Union[str, Path, None] = None,
    return_plot: bool = False,
):
    """(bokeh) Visualize cells in ROI

    Args:
        adata: anndata object
        query: a dict use to select which ROI to display,
            like {"Patients": "Patient 1", "ROI": "ROI3"}, "Patient" and "ROI" are keys in anndata.obs
        geom: "shape" or "point"
        selected_types: select which types to show, others will be muted in grey
        type_key: the key of cell type in anndata.obs (Default: spatialtis.CONFIG.CELL_TYPE_KEY)
        shape_key: the key of cell shape in anndata.obs (Default: spatialtis.CONFIG.SHAPE_KEY)
        centroid_key: the key of cell centroid in anndata.obs (Default: spatialtis.CONFIG.CENTROID_KEY)
        size: size of plot in pixels
        title: title of plot
        palette: config the color, array of color in hex, or
            array of `names of palettes <https://docs.bokeh.org/en/latest/docs/reference/palettes.html>`_
        display: whether to display the plot
        save: path to save your plot
        return_plot: whether to return the plot instance

    """
    if type_key is None:
        type_key = CONFIG.CELL_TYPE_KEY
    if shape_key is None:
        shape_key = CONFIG.SHAPE_KEY
    if centroid_key is None:
        centroid_key = CONFIG.CENTROID_KEY
    if geom not in ["shape", "point"]:
        raise ValueError("Available value for `geom` are 'shape' and 'point'.")

    if geom == "shape":
        if shape_key not in adata.obs.keys():
            geom = "point"
            warnings.warn("Shape key not exist, try to resolve cell as point")

    if geom == "point":
        if centroid_key not in adata.obs.keys():
            raise KeyError("Centroid key not exist")

    df = adata.obs.query("&".join([f"({k}=='{v}')" for k, v in query.items()])).copy()

    if selected_types is not None:
        new_types = []
        for i in df[type_key]:
            if i in selected_types:
                new_types.append(i)
            else:
                new_types.append("other")
        df.loc[:, [type_key]] = new_types

    groups = df.groupby(type_key)

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
        x = [[c[0] for c in eval(cell)] for cell in data[shape_key]]
        y = [[c[1] for c in eval(cell)] for cell in data[shape_key]]
        plot_data = dict(x=x, y=y, name=[name for _ in range(len(x))])
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

    def add_circle(name, fill_color=None, fill_alpha=None):
        cent = [eval(cell) for cell in data[centroid_key]]
        x = [c[0] for c in cent]
        y = [c[1] for c in cent]
        plot_data = dict(x=x, y=y, name=[name for _ in range(len(x))])

        b = p.circle(
            "x",
            "y",
            source=plot_data,
            fill_color=fill_color,
            fill_alpha=fill_alpha,
            line_color="white",
            line_width=0.5,
            size=5,
        )
        if name not in legends_name:
            legends_name.append(name)
            legends.append(LegendItem(label=name, renderers=[b]))

    if selected_types is None:
        for color, (n, data) in zip(colors, groups):
            if geom == "shape":
                add_patches(n, fill_color=color, fill_alpha=0.8)
            else:
                add_circle(n, fill_color=color, fill_alpha=0.8)
    else:
        if geom == "shape":
            for color, (n, data) in zip(colors, groups):
                if n in selected_types:
                    add_patches(n, fill_color=color, fill_alpha=0.8)
                else:
                    add_patches(n, fill_color="grey", fill_alpha=0.5)
        else:
            for color, (n, data) in zip(colors, groups):
                if n in selected_types:
                    add_circle(n, fill_color=color, fill_alpha=0.8)
                else:
                    add_circle(n, fill_color="grey", fill_alpha=0.8)

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
    p.match_aspect = True

    # save something
    if save is not None:
        save_bokeh(p, save)

    # solve env here
    if display is None:
        if CONFIG.WORKING_ENV is None:
            display = False
        else:
            display = True
    if display:
        show(p)

    # it will return a bokeh plot instance, allow user to do some modification
    if return_plot:
        return p
