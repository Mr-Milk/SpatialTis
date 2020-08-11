import warnings
from pathlib import Path
from typing import Mapping, Optional, Sequence, Union

import pyecharts.options as opts
from anndata import AnnData
from pyecharts.charts import Scatter

from spatialtis import CONFIG

from .palette import get_linear_colors


def cell_map_echarts(
    adata: AnnData,
    query: Mapping,
    selected_types: Optional[Sequence] = None,
    type_key: Optional[str] = None,
    centroid_key: Optional[str] = None,
    size: Optional[Sequence[int]] = None,
    title: Optional[str] = None,
    palette: Union[Sequence[str], str, None] = None,
    display: bool = True,
    save: Union[str, Path, None] = None,
    return_plot: bool = False,
):
    """(bokeh) Visualize cells in ROI

    Args:
        adata: anndata object
        query: a dict use to select which ROI to display,
            like {"Patients": "Patient 1", "ROI": "ROI3"}, "Patient" and "ROI" are keys in anndata.obs
        selected_types: select which types to show, others will be muted in grey
        type_key: the key of cell type in anndata.obs (Default: spatialtis.CONFIG.CELL_TYPE_KEY)
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
    if centroid_key is None:
        centroid_key = CONFIG.CENTROID_KEY

    df = adata.obs.query("&".join([f"({k}=='{v}')" for k, v in query.items()]))
    new_types = []
    for i in df[type_key]:
        if i in selected_types:
            new_types.append(i)
        else:
            new_types.append("other")
    df.loc[:, [type_key]] = new_types
    groups = df.groupby(type_key)
    coord = df[centroid_key]

    xdata = [eval(p)[0] for p in coord]

    ydata = dict()
    for t, g in groups:
        ydata[t] = [eval(p)[1] for p in g[centroid_key]]

    initopt_config = dict(width=f"{size[0]}px", height=f"{size[1]}px",)

    a = Scatter(init_opts=opts.InitOpts(**initopt_config))

    a.add_xaxis(xdata)

    for t, data in ydata.items():
        a.add_yaxis(
            str(t),
            data,
            symbol_size=3,
            symbol="rectangle",
            label_opts=opts.LabelOpts(is_show=False),
        )

    a.set_global_opts(
        xaxis_opts=opts.AxisOpts(
            type_="value", is_show=True, is_scale=True, min_interval=1
        ),
        yaxis_opts=opts.AxisOpts(
            type_="value", is_show=True, is_scale=True, min_interval=1
        ),
        title_opts=opts.TitleOpts(),
        tooltip_opts=opts.TooltipOpts(is_show=False),
        toolbox_opts=opts.ToolboxOpts(
            feature={"saveAsImage": {"title": "save", "pixelRatio": 5,},},
        ),
    )

    return a.render_notebook()
