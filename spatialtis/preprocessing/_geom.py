import alphashape
import numpy as np
from shapely.geometry import MultiPoint, Polygon


def geom_cells(cells, method="convex", alpha=0):
    """
    return cell approximate borders, to save computation power

    Args:
        cells: points to describe cells
        method: "convex" using shapely's convex_hull, "concave" using alphashape
        alpha: alpha parameter

    """
    if method == "convex":
        polycells = [MultiPoint(cell).convex_hull for cell in cells]
    elif method == "concave":
        polycells = [alphashape.alphashape(cell, alpha) for cell in cells]
    else:
        raise ValueError("Polygonize options are 'convex' or 'concave'")

    area = []
    borders = []
    centroids = []
    eccentricities = []
    deleted_cells = []
    for i, cell in enumerate(polycells):
        # in case some cell has only one point, or the points become a straight line
        if not isinstance(cell, Polygon):
            deleted_cells.append(i)
        else:
            # cell area
            area.append(cell.area)

            # cell border
            x_border, y_border = cell.exterior.xy
            # [0:-1] to delete the last element, because the start and end of border is the same
            borders.append(tuple((zip(list(x_border), list(y_border))))[0:-1])

            # cell center
            x_ct, y_ct = cell.centroid.xy
            centroids.append((x_ct[0], y_ct[0],))

            # cell eccentricity
            bbox = cell.bounds
            a = np.abs(bbox[0] - bbox[2])
            b = np.abs(bbox[1] - bbox[3])
            ecct = (a / b) if a < b else (b / a)
            eccentricities.append(ecct)

    cells = list(np.delete(np.asarray(cells), deleted_cells))

    return cells, [area, borders, centroids, eccentricities]
