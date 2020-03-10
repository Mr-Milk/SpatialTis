from typing import Union

import numpy as np
from shapely.geometry import box

from ._neighbors import Neighbors


class NeighborsNotFoundError(Exception):
    pass


def check_neighbors(n: Neighbors):
    if not n.neighborsbuilt:
        raise NeighborsNotFoundError(
            "Please run .find_neighbors() before further analysis."
        )
    if n.unitypes is None:
        raise ValueError("The types are not specific.")


def quad_tessellation(
    rect: tuple, quad: Union[tuple, str] = "auto", grid_size: float = 150,
):
    if quad == "auto":
        nx = int((rect[2] - rect[0]) / grid_size)
        ny = int((rect[3] - rect[1]) / grid_size)
        quad = (nx, ny)

    xrange = np.linspace(rect[0], rect[2], quad[0] + 1)
    yrange = np.linspace(rect[1], rect[3], quad[1] + 1)

    rects = []
    for i, x in enumerate(xrange[0:-1]):
        x1 = x
        x2 = xrange[i + 1]
        lower = [(x1, y) for y in yrange[0:-1]]
        higher = [(x2, y) for y in yrange[1::]]
        for p in tuple(zip(lower, higher)):
            rects.append(box(*tuple(np.asarray(p).ravel())))
    return rects, quad
