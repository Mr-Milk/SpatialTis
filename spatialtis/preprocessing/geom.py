import numpy as np
from shapely.geometry import MultiPoint


def geom_cells(cells):
    """
    return cell approximate borders, to save computation power
    """
    polycells = [MultiPoint(cell).convex_hull for cell in cells]

    area = []
    borders = []
    centroids = []
    eccentricities = []
    for cell in polycells:
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

    # print(len(area), len(borders), len(centroids), len(eccentricities))

    return [area, borders, centroids, eccentricities]
