from collections import Counter
from typing import Optional, Union, Sequence, Mapping
from pathlib import Path

import pandas as pd

from pyecharts import options as opts
from pyecharts.charts import Pie

from ._save import save_pyecharts
from .palette import get_linear_colors


def get_grid(counts, min_ncol=8):
    if counts < min_ncol:
        return 1, counts
    else:
        nrow = counts // min_ncol
        exceed = counts % min_ncol
        if exceed:
            return nrow + 1, min_ncol
        else:
            return nrow, min_ncol


def grouped_pie(df: pd.DataFrame,
                mapper: Mapping,
                order: Optional[Sequence] = None,
                selected_types: Optional[Sequence] = None,
                pie_size: float = 0.8,
                round_size: float = 0.7,
                size: Sequence = (900, 500),
                renderer: str = 'canvas',
                theme: str = 'white',
                palette: Optional[Sequence] = None,
                display: bool = True,
                return_plot: bool = False,
                title: Optional[str] = None,
                save: Union[str, Path, None] = None,
                ):
    """

    Args:
        df:
        mapper:
        order:
        selected_types:
        pie_size:
        round_size:
        size:
        renderer:
        theme:
        palette:
        display:
        return_plot:
        title:
        save:

    Returns:

    """
    if palette is not None:
        palette = get_linear_colors(palette)

    if selected_types is not None:
        df = df[selected_types]

    c = Pie(
        init_opts=opts.InitOpts(width=f"{size[0]}px",
                                height=f"{size[1]}px",
                                renderer=renderer,
                                theme=theme,
                                )
    )

    width = int(''.join([i for i in c.width if i.isdigit()]))
    height = int(''.join([i for i in c.height if i.isdigit()]))

    (nrow, ncol) = get_grid(df.shape[1])

    grid_width = width / ncol
    grid_height = height*0.95 / nrow

    step_w = grid_width / width * 100
    start_w = step_w / 2

    step_h = grid_height / height * 100
    start_h = step_h / 2 + 5

    center_w = [str(start_w + step_w * i) + '%' for i in range(ncol)]
    center_h = [str(start_h + step_h * i) + '%' for i in range(nrow)]

    centers = []
    for w in center_w:
        for h in center_h:
            centers.append([w, h])

    if grid_width <= grid_height:
        outer_r = grid_width / 2 * pie_size
    else:
        outer_r = grid_height / 2 * pie_size
    inner_r = outer_r * round_size
    r = [inner_r, outer_r]

    unitypes = pd.unique(df.values.ravel())

    for i, (label, data) in enumerate(df.items()):
        c_data = Counter(data)
        for k in unitypes:
            if k not in c_data.keys():
                c_data[k] = 0
        if order is None:
            order = list(c_data.keys())

        c.add(
            label,
            [(mapper[order[i]], c_data[order[i]]) for i in range(len(order))],
            center=centers[i],
            radius=r,
            rosetype="radius",
            label_opts=opts.LabelOpts(formatter="{a}", position="center"),
        )

    c.set_global_opts(
        title_opts=opts.TitleOpts(title=title),
        # visualmap_opts=opts.VisualMapOpts(range_color=palette),
        toolbox_opts=opts.ToolboxOpts(feature={"saveAsImage": {"title": "save", "pixelRatio": 5, }, }, ),
    )

    if save is not None:
        save_pyecharts(c, save)

    if display:
        c.load_javascript()
        return c.render_notebook()

    if return_plot:
        return c
