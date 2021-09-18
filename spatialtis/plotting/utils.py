from itertools import cycle
from typing import Union, List

from matplotlib import cm
from matplotlib.colors import Colormap
from matplotlib.colors import to_hex

COLOR_POOL = cycle(['#1f77b4ff', '#aec7e8ff', '#ff7f0eff', '#ffbb78ff',
                    '#2ca02cff', '#98df8aff', '#d62728ff', '#ff9896ff',
                    '#9467bdff', '#c5b0d5ff', '#8c564bff', '#c49c94ff',
                    '#e377c2ff', '#f7b6d2ff', '#bcbd22ff', '#dbdb8dff',
                    '#17becfff', '#9edae5ff', '#e41a1cff', '#377eb8ff',
                    '#4daf4aff', '#984ea3ff', '#ff7f00ff', '#ffff33ff',
                    '#a65628ff', '#f781bfff'])
