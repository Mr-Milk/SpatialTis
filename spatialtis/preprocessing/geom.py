from pathlib import Path
from typing import Optional, Sequence, Union

import numpy as np
from shapely.geometry import MultiPoint


def mask2cells(
    mask_img: Union[Path, str],
    bg: Optional[int] = None,
    polygonize: str = "convex",
    alpha: Optional[float] = None,
):
    try:
        from skimage.io import imread
        from skimage.measure import label, regionprops
    except ImportError:
        raise ImportError(
            "Using preprocessing moduele needs scikit-image. "
            "Try pip install scikit-image."
        )

    area = []
    borders = []
    centroids = []
    eccentricities = []

    cells = []

    mask = imread(mask_img)
    label_mask = label(mask, background=bg)

    for cell in regionprops(label_mask):
        border = cell_border(cell.coords, polygonize=polygonize, alpha=alpha)
        if border is not None:
            x_border, y_border = border.exterior.xy

            area.append(cell.area)
            # [0:-1] to delete the last element, because the start and end of border is the same
            borders.append(tuple((zip(list(x_border), list(y_border))))[0:-1])
            centroids.append(cell.centroid)
            eccentricities.append(cell.eccentricity)
            cells.append(cell.coords)

    return cells, [area, borders, centroids, eccentricities]


def cell_border(cell, polygonize="convex", alpha=None):
    # if cell has less than 3 points, then it's meaningless
    if len(cell) < 3:
        return None
    if polygonize == "convex":
        cell = MultiPoint(cell).convex_hull
    elif polygonize == "concave":
        try:
            import alphashape
        except ImportError:
            raise ImportError(
                "You need alphashape installed, " "try `pip install alphashape`."
            )
        if alpha is None:
            alpha = alphashape.optimizealpha(cell)
        try:
            cell = alphashape.alphashape(cell, alpha)
        except Exception:
            return None
    else:
        raise ValueError("Polygonize options are 'convex' or 'concave'")
    if not hasattr(cell, "exterior"):
        return None
    return cell


def get_cell_exp_stack(
    stack: Sequence, cells: Sequence, method: str = "mean"
) -> Sequence:
    """
    Parameters
    channel: list or ndarray, matrix info of the channel
    cells: list or ndarray, each element contains points for each cell
    method: str, ('mean' / 'median' / 'sum' / ...)  any numpy method to compute expression level of single cell
    """
    cells_density = list()
    for cell in cells:
        cell_pixels = stack[:, [i[0] for i in cell], [i[1] for i in cell]]
        exec(f"cells_density.append([np.{method}(density) for density in cell_pixels])")
    # each secondary array as a channel
    return cells_density
