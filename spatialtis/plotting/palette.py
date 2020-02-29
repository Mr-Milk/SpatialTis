import logging
from itertools import cycle

import bokeh.palettes as pl
import numpy as np
import seaborn as sns
from colour import Color


def colorcycle(*palette):
    new_palette = list()
    palette_keys = pl.all_palettes.keys()
    for p in palette:
        if p in palette_keys:
            pcolors: dict = pl.all_palettes[p]
            max_color = pcolors[max(pcolors.keys())]
            for c in max_color:
                if c not in new_palette:
                    new_palette.append(c)
        else:
            try:
                c = Color(p).hex_l
                if c not in new_palette:
                    new_palette.append(c)
            except ValueError:
                raise ValueError(f"'{p}' is not a palette name nor a color")

    new_palette = sorted(set(new_palette), key=new_palette.index)

    return cycle(new_palette), new_palette


def get_colors(n: int, palette):
    cycler, _ = colorcycle(*palette)
    return [next(cycler) for i in range(0, n)]


def get_linear_colors(palette):
    cycler, colors = colorcycle(*palette)
    max_colors = len(colors)
    color_index = np.linspace(0, max_colors - 1, max_colors)
    return [colors[int(i)] for i in color_index]


def view_palette(palette):
    return sns.palplot(palette)
