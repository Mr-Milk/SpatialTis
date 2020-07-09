from typing import Mapping, Optional, Sequence

import pyecharts.options as opts
from anndata import AnnData
from pyecharts.charts import Bar3D, Scatter, Tab

from ..config import CONFIG
from .palette import get_linear_colors


def expression_map(
    adata: AnnData,
    query: Mapping,
    marker_key: Optional[str] = None,
    centroid_key: Optional[str] = None,
    order: Optional[Sequence] = None,
    use: str = "bar3d",  # 'bar3d', 'scatter'
    renderer: str = "canvas",
    axis_size: tuple = (100, 100, 80),
    size: tuple = (800, 500),
    palette: Optional[Sequence] = None,
    display: bool = True,
    # save: Union[str, Path, None] = None, # save multi plots is not allowed in pyecharts
    return_plot: bool = False,
):
    """(pyecharts) Visualize marker expression in ROI

    Issue: There are problems in saving this plot even in .html
        For now please use save bottom.

    Args:
        adata: anndata object
        query: a dict use to select which ROI to display,
            like {"Patients": "Patient 1", "ROI": "ROI3"}, "Patient" and "ROI" are keys in anndata.obs
        marker_key: the key of marker in anndata.var (Default: spatialtis.CONFIG.MARKER_KEY)
        centroid_key: the key of cell centroid in anndata.obs (Default: spatialtis.CONFIG.CENTROID_KEY)
        order: array of marker name, display as order
        use: "bar3d" or "scatter"
        renderer: "canvas" or "svg"
        axis_size: the length of x,y,z axis
        size: size of plot in pixels
        palette: config the color, array of color in hex, or
            array of `names of palettes <https://docs.bokeh.org/en/latest/docs/reference/palettes.html>`_
        display: whether to display the plot
        return_plot: whether to return the plot instance

    """
    if marker_key is None:
        marker_key = CONFIG.MARKER_KEY
    if centroid_key is None:
        centroid_key = CONFIG.CENTROID_KEY
    if use not in ["bar3d", "scatter"]:
        raise ValueError(
            "No such plot method, available options are 'bar3d' and 'scatter'."
        )

    if palette is None:
        palette = ["RdYlBu"]
    default_color = get_linear_colors(palette)

    gene_names = list(adata.var[marker_key])
    data = adata.obs.query("&".join([f"({k}=='{v}')" for k, v in query.items()]))
    coord = data[centroid_key]
    if order is not None:
        exp_index = [gene_names.index(i) for i in order]
        gene_names = order
        exp = adata[data.index].X.T[exp_index]
    else:
        # if not specific, it will take the first 3
        exp = adata[data.index].X.T[0:3]

    t = Tab()
    for iexp, gene_name in zip(exp, gene_names):
        zdata = []

        for e, c in zip(iexp, coord):
            p = eval(c)
            zdata.append([p[0], p[1], float(e)])

        zrange = sorted(zdata, key=lambda k: k[2])
        initopt_config = dict(
            width=f"{size[0]}px", height=f"{size[1]}px", renderer=renderer,
        )
        if use == "bar3d":
            a = Bar3D(init_opts=opts.InitOpts(**initopt_config))

            a.add(
                series_name="",
                shading="color",
                data=zdata,
                xaxis3d_opts=opts.Axis3DOpts(type_="value"),
                yaxis3d_opts=opts.Axis3DOpts(type_="value"),
                zaxis3d_opts=opts.Axis3DOpts(type_="value"),
                grid3d_opts=opts.Grid3DOpts(
                    width=axis_size[1], height=axis_size[2], depth=axis_size[0]
                ),
            ).set_global_opts(
                title_opts=opts.TitleOpts(gene_name),
                visualmap_opts=opts.VisualMapOpts(
                    dimension=2,
                    max_=zrange[-1][2],
                    min_=zrange[0][2],
                    range_color=default_color,
                ),
                tooltip_opts=opts.TooltipOpts(is_show=False),
                toolbox_opts=opts.ToolboxOpts(
                    feature={"saveAsImage": {"title": "save", "pixelRatio": 5,},},
                ),
            )
        else:
            a = Scatter(init_opts=opts.InitOpts(**initopt_config))

            a.add_xaxis([i[0] for i in zdata]).add_yaxis(
                "",
                [i[1::] for i in zdata],
                symbol_size=3,
                label_opts=opts.LabelOpts(is_show=False),
            ).set_global_opts(
                xaxis_opts=opts.AxisOpts(
                    type_="value", is_show=True, is_scale=True, min_interval=1
                ),
                yaxis_opts=opts.AxisOpts(
                    type_="value", is_show=True, is_scale=True, min_interval=1
                ),
                title_opts=opts.TitleOpts(gene_name),
                tooltip_opts=opts.TooltipOpts(is_show=False),
                toolbox_opts=opts.ToolboxOpts(
                    feature={"saveAsImage": {"title": "save", "pixelRatio": 5,},},
                ),
                visualmap_opts=opts.VisualMapOpts(
                    type_="color",
                    min_=zrange[0][2],
                    max_=zrange[-1][2],
                    range_color=default_color,
                    dimension=2,
                ),
            )

        t.add(a, gene_name)

    '''
    if save is not None:
        # nested tab can only save in html
        p = Path(save)
        if p.suffix[1:] != 'html':
            p += '.html'
        # save_path = f"""{'/'.join(p.parts[:-1])}/{p.stem}.html"""
        t.render(p)
    '''

    if display:
        t.load_javascript()
        return t.render_notebook()

    if return_plot:
        return t
