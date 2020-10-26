from collections import OrderedDict
from typing import Any

from shapely.geometry import MultiPoint

from .neighbors import Neighbors


class NeighborsNotFoundError(Exception):
    pass


def check_neighbors(n: Any):
    if not isinstance(n, Neighbors):
        raise TypeError("A spatialtis.Neighbors subject is needed.")

    if not n.neighborsbuilt:
        raise NeighborsNotFoundError(
            "Please run .find_neighbors() before further analysis."
        )
    if n.unitypes is None:
        raise ValueError("The types are not specific.")


# modify from PySAL pointpats
class QuadStats:
    def __init__(self, points, nx=None, ny=None, grid_size=None):

        self.points = points
        self.bbox = MultiPoint(points).bounds
        self.width = self.bbox[2] - self.bbox[0]
        self.height = self.bbox[3] - self.bbox[1]

        if (nx is None) & (ny is None):
            self.nx = int(self.width // grid_size)
            self.ny = int(self.height // grid_size)
        else:
            self.nx = nx
            self.ny = ny

        if (self.nx != 0) & (self.ny != 0):
            self.w_x = self.width / self.nx
            self.h_y = self.height / self.ny

        self.cells_grid_id = []

    def grid_counts(self):
        dict_id_count = OrderedDict()
        for i in range(self.ny):
            for j in range(self.nx):
                dict_id_count[j + i * self.nx] = 0
        for point in self.points:
            index_x = int((point[0] - self.bbox[0]) // self.w_x)
            index_y = int((point[1] - self.bbox[1]) // self.h_y)
            if index_x == self.nx:
                index_x -= 1
            if index_y == self.ny:
                index_y -= 1
            id_ = index_y * self.nx + index_x
            self.cells_grid_id.append(id_)
            dict_id_count[id_] += 1

        return dict_id_count
