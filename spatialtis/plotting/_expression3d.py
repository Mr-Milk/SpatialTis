from typing import Union, Optional, Sequence, Mapping
from pathlib import Path

import pandas as pd
from anndata import AnnData

import pyecharts.options as opts
from pyecharts.charts import Bar3D, Scatter3D, Tab

from ..config import CONFIG
from .palette import get_linear_colors


def expression_map(
                    adata: AnnData,
                    query: Mapping,
                    marker_col: Optional[str] = None,
                    centroid_col: Optional[str] = None,
                    order: Optional[Sequence] = None,
                    method: str = 'bar3d',  # 'bar3d', 'scatter3d'
                    renderer: str = 'canvas',
                    axis_size: Sequence = (100, 100, 80),
                    size: Sequence = (800, 500),
                    palette: Optional[Sequence] = None,
                    display: bool = True,
                    # save: Union[str, Path, None] = None, # save multi plots is not allowed in pyecharts
                    return_plot: bool = False,
                    ):
    if marker_col is None:
        marker_col = CONFIG.MARKER_COL
    if centroid_col is None:
        centroid_col = CONFIG.CENTROID_COL
    if method not in ['bar3d', 'scatter3d']:
        raise ValueError("No such plot method, available options are 'bar3d' and 'scatter3d'.")

    if palette is None:
        palette = ['RdYlBu']
    default_color = get_linear_colors(palette)

    gene_names = list(adata.var[marker_col])
    data = adata.obs.query("&".join([f"({k}=='{v}')" for k, v in query.items()]))
    coord = data[centroid_col]
    if order is not None:
        exp_index = [gene_names.index(i) for i in order]
        gene_names = order
        exp = adata[data.index].X.T[exp_index]
    else:
        # if not specific, it will take the first 10
        exp = adata[data.index].X.T[0:10]

    t = Tab()
    for x, iexp in enumerate(exp):
        gene_name = gene_names[x]
        zdata = []

        for i, c in enumerate(coord):
            p = eval(c)
            zdata.append([p[0], p[1], iexp[i]])

        zrange = sorted(zdata, key=lambda k: k[2])
        initopt_config = dict(
            width=f"{size[0]}px",
            height=f"{size[1]}px",
            renderer=renderer,
        )
        if method == 'bar3d':
            a = Bar3D(init_opts=opts.InitOpts(**initopt_config))
        else:
            a = Scatter3D(init_opts=opts.InitOpts(**initopt_config))

        a.add(
            series_name="",
            shading="color",
            data=zdata,
            xaxis3d_opts=opts.Axis3DOpts(type_="value"),
            yaxis3d_opts=opts.Axis3DOpts(type_="value"),
            grid3d_opts=opts.Grid3DOpts(width=axis_size[1], height=axis_size[2], depth=axis_size[0]),
        ).set_global_opts(
            visualmap_opts=opts.VisualMapOpts(
                dimension=2,
                max_=zrange[-1][2],
                min_=zrange[0][2],
                range_color=default_color,
            ),
            tooltip_opts=opts.TooltipOpts(is_show=False),
            toolbox_opts=opts.ToolboxOpts(feature={
                "saveAsImage": {"title": "save", "pixelRatio": 5, },
            }, )
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
        # t.load_javascript()
        return t.render_notebook()

    if return_plot:
        return t
