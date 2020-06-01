from pathlib import Path
from typing import Optional, Sequence, Union

from pyecharts import options as opts
from pyecharts.charts import Sankey

from ._save import save_pyecharts


def sankey(
    nodes: Sequence,
    nodes_colors: Sequence,
    links: Sequence,
    renderer: str = "canvas",
    theme: str = "white",
    size: Sequence = (900, 500),
    title: Optional[str] = None,
    display: bool = True,
    save: Union[str, Path, None] = None,
    return_plot: bool = False,
):
    """

    Args:
        nodes:
        nodes_colors:
        links:
        renderer:
        theme:
        size:
        title:
        display:
        save:
        return_plot:

    Returns:

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
            None,
            draw_nodes,
            draw_links,
            linestyle_opt=opts.LineStyleOpts(color="source", opacity=0.6, curve=0.5),
            label_opts=opts.LabelOpts(position="right"),
        )
        .set_global_opts(
            title_opts=opts.TitleOpts(title=title),
            toolbox_opts=opts.ToolboxOpts(
                feature={"saveAsImage": {"title": "save", "pixelRatio": 5,},},
            ),
        )
    )

    if save is not None:
        save_pyecharts(c, save)

    if display:
        c.load_javascript()
        return c.render_notebook()

    if return_plot:
        return c
