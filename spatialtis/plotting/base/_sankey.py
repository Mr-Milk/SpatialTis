from pathlib import Path
from typing import Optional, Sequence, Union

from pyecharts import options as opts
from pyecharts.charts import Sankey

from spatialtis.config import CONFIG
from spatialtis.utils import reuse_docstring

from .save import save_pyecharts


@reuse_docstring()
def sankey(
    nodes: Sequence,
    nodes_colors: Sequence,
    links: Sequence,
    renderer: str = "canvas",
    theme: str = "white",
    size: Sequence = (900, 500),
    title: Optional[str] = None,
    display: Optional[bool] = None,
    save: Union[str, Path, None] = None,
    return_plot: bool = False,
):
    """(pyecharts) Sankey plot

    Args:
        nodes: All the nodes
        nodes_colors: The color for each node
        links: Connections between nodes
        renderer: {renderer}
        theme: {theme}
        size: {size}
        title: {title}
        display: {display}
        save: {save}
        return_plot: {return_plot}

    """
    draw_nodes = []
    draw_links = []
    for n, c in zip(nodes, nodes_colors):
        draw_nodes.append(
            {"name": n, "itemStyle": {"normal": {"color": c, "borderColor": "#aaa"}}}
        )

    for (source, target, value) in links:
        draw_links.append({"source": source, "target": target, "value": value})

    c = (
        Sankey(
            init_opts=opts.InitOpts(
                width=f"{size[0]}px",
                height=f"{size[1]}px",
                renderer=renderer,
                theme=theme,
            )
        )
        .add(
            "",
            draw_nodes,
            draw_links,
            linestyle_opt=opts.LineStyleOpts(color="source", opacity=0.6, curve=0.5),
            label_opts=opts.LabelOpts(position="right"),
        )
        .set_global_opts(
            title_opts=opts.TitleOpts(title=title),
            toolbox_opts=opts.ToolboxOpts(
                feature={"saveAsImage": {"title": "save", "pixelRatio": 5}},
            ),
        )
    )

    if save is not None:
        save_pyecharts(c, save)

    if display is None:
        if CONFIG.WORKING_ENV is None:
            display = False
        else:
            display = True
    if display:
        c.load_javascript()
        return c.render_notebook()

    if return_plot:
        return c
