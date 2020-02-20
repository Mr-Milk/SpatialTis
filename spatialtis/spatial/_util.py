from typing import Sequence

from ._neighbors import Neighbors


class NeighborsNotFoundError(Exception):
    pass


def check_neighbors(n: Neighbors):
    if not n.treebuilt:
        raise NeighborsNotFoundError("Please run .find_neighbors() before further analysis.")
    if n.types is None:
        raise ValueError("The types are not specific.")
