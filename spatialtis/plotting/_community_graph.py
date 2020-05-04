from collections import Counter
from typing import Optional, Union, Sequence
from pathlib import Path

import pandas as pd

from pyecharts import options as opts
from pyecharts.charts import Graph

from ._save import save_pyecharts
from .palette import get_linear_colors, get_colors


def graph_plot(df: pd.DataFrame,
               node_col: Optional[str] = None,  # centroid_col
               node_category_col: Optional[str] = None,  # cell type / communities
               node_info_col: Optional[str] = None,
               edge_col: Optional[str] = None,  # neighbors relationship
               edge_category_col: Optional[str] = None,
               edge_info_col: Optional[str] = None,
               node_size: Union[float, int] = 3,
               edge_size: Union[float, int] = 0.1,
               size: Sequence = (800, 800),
               renderer: str = 'canvas',
               theme: str = 'white',
               palette: Optional[Sequence] = None,
               display: bool = True,
               return_plot: bool = False,
               title: Optional[str] = None,
               save: Union[str, Path, None] = None,
               ):

    if palette is not None:
        palette = get_linear_colors(["Set3"])

    nodes_data = []
    edges_data = []
    categories = []

    cols = list(df.columns)
    ixy = cols.index(node_col)
    iedge = cols.index(edge_col)
    if node_category_col is not None:
        inode_category = cols.index(node_category_col)
    if node_info_col is not None:
        inode_info = cols.index(node_info_col)
    if edge_category_col is not None:
        iedge_category = cols.index(edge_category_col)
        edge_categories = df[edge_category_col]
        edge_types = pd.unique(edge_categories)
        edges_colors = dict(zip(edge_types, get_colors(len(edge_types), ["Set3", "Spectral"])))
    if edge_info_col is not None:
        iedge_info = cols.index(edge_info_col)

    for i, (_, c) in enumerate(df.iterrows()):
        xy = eval(c[ixy])
        node_config = dict(
            name=str(i),
            x=xy[1],
            y=xy[0],
            label_opts=opts.LabelOpts(is_show=False),
            symbol_size=node_size,
        )
        if node_category_col is not None:
            category = str(c[inode_category])
            node_config['category'] = category
            categories.append(opts.GraphCategory(name=category))
        if node_info_col is not None:
            node_config['value'] = str(c[inode_info])

        nodes_data.append(opts.GraphNode(**node_config))

        for n in c[iedge]:
            if edge_category_col is not None:
                source_category = edge_categories[n]
                target_category = edge_categories[i]
                if source_category == target_category:
                    edges_data.append(opts.GraphLink(source=str(n), target=str(i), linestyle_opts=opts.LineStyleOpts(
                        width=edge_size, color=edges_colors[source_category],
                    )))
            else:
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
        # visualmap_opts=opts.VisualMapOpts(range_color=palette),
        legend_opts=opts.LegendOpts(type_="scroll", orient="vertical", pos_left="2%", pos_top="20%"),
        toolbox_opts=opts.ToolboxOpts(feature={
            "saveAsImage": {"title": "save", "pixelRatio": 5, },
            "brush": {},
            "restore": {},
        }, ),
    )

    if save is not None:
        save_pyecharts(g, save)

    if display:
        g.load_javascript()
        return g.render_notebook()

    if return_plot:
        return g
