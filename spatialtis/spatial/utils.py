from ast import literal_eval
from collections import OrderedDict

import numpy as np


# modify from PySAL pointpats
class QuadStats:
    def __init__(self, points, bbox, nx=None, ny=None, grid_size=None):

        self.points = points
        self.bbox = bbox
        minx, miny, maxx, maxy = self.bbox
        self.width = maxx - minx
        self.height = maxy - miny

        if (nx is None) & (ny is None):
            self.nx = int(self.width // grid_size)
            self.ny = int(self.height // grid_size)
            self.nx = self.nx if self.nx != 0 else 1
            self.ny = self.ny if self.ny != 0 else 1
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


def job_cutter(total, nums):
    if total > nums:
        order = int(total / nums) * nums
        if order != total:
            k = [i for i in range(0, order, int(order / nums))] + [order, total]
        else:
            k = [i for i in range(0, order, int(order / nums))] + [total]
    else:
        k = [0, total]
    lk = len(k)
    return [(k[i], k[i + 1]) for i in range(lk) if i < (lk - 1)]


def normalize(arr: np.ndarray):
    min_max = arr.amax() - arr.amin()
    if min_max != 0:
        return (arr - arr.amin()) / min_max
    else:
        return np.ones(arr.shape)


def get_eval(df, key, need_eval):
    if need_eval:
        return [literal_eval(c) for c in df[key]]
    else:
        return [c for c in df[key]]
