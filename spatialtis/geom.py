import numpy as np
from shapely.geometry import Polygon, MultiPoint, Point
from shapely.strtree import STRtree
from itertools import combinations
from scipy.spatial import cKDTree
from collections import Counter

def border(cell):
    return MultiPoint(cell).convex_hull


def polygonize_cells(cells):
    '''
    return cell approximate borders, to save computation power
    '''
    border = [MultiPoint(cell).convex_hull.exterior.xy for cell in cells]
    # the start point and the end point of a convex_hull polygon are the same, just get rid of the end point
    return [tuple(zip(list(x),list(y)))[0:-1] for (x,y) in border]

def relationships_kdtree(polycells, expand_length=5):
    '''
    using kd-tree
    '''
    for ic, c in enumerate(polycells):
        c.ix = ic
    count = len(polycells)
    relationships = dict(zip(range(0, count), [[] for i in range(0, count)]))
    centroids = [(c.x, c.y) for c in [cell.centroid for cell in polycells]]
    diameter = np.median([np.sqrt(cell.area) for cell in polycells])
    tree = cKDTree(centroids)
    pairs = []
    for ix, point in enumerate(centroids):
        result = tree.query_ball_point(point, 2*diameter+expand_length)
        for i in result:
            if ix != i:
                pairs.append((ix, i))
    pairs = Counter(frozenset(ele) for ele in pairs).keys()
    
    # TODO: scale the objects, and then test intersects
    for (p1, p2) in pairs:
        if polycells[p1].intersects(polycells[p2]):
            relationships[p1].append(p2)
            relationships[p2].append(p1)

    return relationships

def relationships_Rtree(polycells, scale_factor=1.2):
    '''
    Using shapely in-built STR-packed R-tree
    '''
    count = len(polycells)
    relationships = dict(zip(range(0, count), [[] for i in range(0, count)]))
    # create STR R*tree
    tree = STRtree(polycells)
    for ic, cell in enumerate(polycells):
        # expand cell region
        scaled_cell = scale(cell, xfact=scale_factor, yfact=scale_factor)
        result = tree.query(cell)
        for n in result:
            if n.ix != ic:
                relationships[ic].append(n.ix)
    
    return relationships

def distribution_score():
    pass

def heterogenity_score():
    pass

