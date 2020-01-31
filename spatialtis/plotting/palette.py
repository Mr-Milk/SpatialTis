from itertools import cycle

import bokeh.palettes as pl

import seaborn as sns


def colorcycle(*palette: str):
    new_palette = list()
    for p in palette:
        pcolors: dict = pl.all_palettes[p]
        max_color = pcolors[max(pcolors.keys())]
        new_palette += max_color

    return cycle(new_palette)


def get_colors(cycler: cycle, n: int):
    return [next(cycler) for i in range(0, n)]


def view_palette(palette):
    return sns.palplot(palette)
