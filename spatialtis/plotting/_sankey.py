from pathlib import Path
from typing import Optional, Sequence, Union

from pyecharts import options as opts
from pyecharts.charts import Sankey

from ..config import CONFIG
from ._save import save_pyecharts


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
    """(pyecharts) sankey plot

    Args:
        nodes: all the nodes
        nodes_colors: the color for each node
        links: connections between nodes
        renderer: "canvas" or "svg"
        theme: https://pyecharts.org/#/zh-cn/themes
        size: size of plot in pixels
        title: title of the plot
        display: whether to display the plot
        save: the path to save your plot
        return_plot: whether to return the plot instance

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
