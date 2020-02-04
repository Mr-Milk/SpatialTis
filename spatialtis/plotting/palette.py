from itertools import cycle

import bokeh.palettes as pl

import seaborn as sns

from colour import Color


def colorcycle(*palette):
    new_palette = list()
    palette_keys = pl.all_palettes.keys()
    for p in palette:
        if p in palette_keys:
            pcolors: dict = pl.all_palettes[p]
            max_color = pcolors[max(pcolors.keys())]
            new_palette += max_color
        else:
            try:
                c = Color(p).hex_l
                new_palette.append(c)
            except ValueError:
                raise ValueError(f"'{p}' is not a palette name nor a color")

    return cycle(new_palette)


def get_colors(n: int, *palette):
    cycler = colorcycle(*palette)
    return [next(cycler) for i in range(0, n)]


def view_palette(palette):
    return sns.palplot(palette)
