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
