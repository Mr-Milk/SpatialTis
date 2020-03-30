from collections import Counter
from typing import Optional, Union, Sequence
from pathlib import Path

import pandas as pd

from pyecharts import options as opts
from pyecharts.charts import Graph

from ._save import save_pyecharts
from .palette import get_linear_colors
from ..config import CONFIG


def graph_plot(df: pd.DataFrame,
               type_col: Optional[str] = None,
               neighbors_col: Optional[str] = None,
               centroid_col: Optional[str] = None,
               tooltip_col: Optional[str] = None,
               size: Sequence = (800, 800),
               renderer: str = 'canvas',
               theme: str = 'white',
               palette: Optional[Sequence] = None,
               display: bool = True,
               return_plot: bool = False,
               title: Optional[str] = None,
               save: Union[str, Path, None] = None,
               ):
    if type_col is None:
        type_col = CONFIG.CELL_TYPE_COL
    if neighbors_col is None:
        neighbors_col = CONFIG.NEIGHBORS_COL
    if centroid_col is None:
        centroid_col = CONFIG.CENTROID_COL
    if palette is not None:
        palette = get_linear_colors(palette)

    nodes_data = []
    edges_data = []
    categories = []

    cols = list(df.columns)
    ixy = cols.index(centroid_col)
    icategory = cols.index(type_col)
    ineighbors = cols.index(neighbors_col)

    if tooltip_col is not None:
        iname = cols.index(type_col)
    else:
        iname = None

    for i, (_, c) in enumerate(df.iterrows()):
        xy = eval(c[ixy])
        category = str(c[icategory])
        node_config = dict(
            name=str(i),
            x=xy[1],
            y=xy[0],
            category=category,
            label_opts=opts.LabelOpts(is_show=False),
        )
        if iname is not None:
            node_config['value'] = str(c[iname])

        nodes_data.append(opts.GraphNode(**node_config))
        categories.append(opts.GraphCategory(name=category))

        for n in c[ineighbors]:
            edges_data.append(opts.GraphLink(source=str(n), target=str(i)))

    g = Graph(init_opts=opts.InitOpts(width=f"{size[0]}px",
                                      height=f"{size[1]}px",
                                      renderer=renderer,
                                      theme=theme,
                                      ))
    g.add("",
          nodes_data,
          edges_data,
          categories,
          layout="none",
          is_rotate_label=True,
          edge_label=opts.LabelOpts
          (is_show=False),
          tooltip_opts=opts.TooltipOpts
          (formatter="{c}"),
          ).set_global_opts(
        title_opts=opts.TitleOpts(title=title),
        visualmap_opts=opts.VisualMapOpts(range_color=palette),
        legend_opts=opts.LegendOpts(type_="scroll", orient="vertical", pos_left="2%", pos_top="20%"),
        toolbox_opts=opts.ToolboxOpts(feature={
            "saveAsImage": {"title": "save", "pixelRatio": 5, },
        }, ),
    )

    if save is not None:
        save_pyecharts(g, save)

    if display:
        g.load_javascript()
        return g.render_notebook()

    if return_plot:
        return g
