import numpy as np
import pandas as pd
from anndata import AnnData
from itertools import combinations_with_replacement
from spatialtis_core import multipoints_bbox, multipolygons_area, polygons_area

from spatialtis.abc import AnalysisBase
from spatialtis.utils import col2adata, doc, read_shapes


@doc
def cell_components(
        data: AnnData,
        export_key: str = "cell_components",
        **kwargs,
):
    """Count the proportion of each types of cells in each group.

    Parameters
    ----------
    data : {adata}
    export_key : {export_key}
    **kwargs : {analysis_kwargs}

    """
    ab = AnalysisBase(data, display_name="Cell components", export_key=export_key, **kwargs)
    ab.check_cell_type()
    result = ab.type_counter()
    result.columns.name = 'cell type'
    ab.result = result


@doc
def cell_density(data: AnnData,
                 ratio: float = 1.0,
                 export_key: str = "cell_density",
                 **kwargs):
    """Calculating cell density in each ROI.

    The size of each ROI will be auto-computed as convex hull of all the cells in a ROI.

    Parameters
    ----------
    data : {adata}
    ratio : float, default: 1.0
        The ratio between the unit used in your dataset and real length unit.
        ratio = Dataset unit / real length unit.
        For example, if the resolution of your dataset is 1μm, but you want to use 1mm as unit,
        then you should set the ratio as 0.001, 1 pixel = 0.001mm.
    export_key : {export_key}
    **kwargs : {analysis_kwargs}

    """
    ab = AnalysisBase(data, display_name="Cell density", export_key=export_key, **kwargs)
    ab.check_cell_type()
    result = ab.type_counter()

    area = []
    for roi_name, points in ab.iter_roi(fields=['centroid']):
        area.append(polygons_area(points))

    area = np.asarray(area) * (ratio * ratio)
    result = result.div(area, axis=0)
    result.columns.name = 'cell type'
    ab.result = result


def _bbox_eccentricity(bbox) -> float:
    x = (bbox[2] - bbox[0]) / 2.0
    y = (bbox[3] - bbox[1]) / 2.0
    if x < y:
        x, y = y, x
    return np.sqrt(1.0 - y ** 2 / x ** 2)


@doc
def cell_morphology(data: AnnData,
                    area_key: str = None,
                    eccentricity_key: str = None,
                    **kwargs):
    """Cell morphology variation between different groups.

    This function only works for data with cell shape information.
    The area is calculated using shoelace formula
    The eccentricity is assumed that the cell is close to ellipse,
    the semi-minor and semi-major axis
    is get from the bbox side.

    Parameters
    ----------
    data : {adata}
    area_key : str
        The `obs` key to store cell area value.
    eccentricity_key : str
        The `obs` key to store cell eccentricity.
    **kwargs : {analysis_kwargs}

    """
    ab = AnalysisBase(data, display_name="Cell morphology", **kwargs)
    shapes = read_shapes(data.obs, ab.shape_key)
    areas = multipolygons_area(shapes)
    eccentricity = [_bbox_eccentricity(bbox) for bbox in multipoints_bbox(shapes)]
    area_key = ab.area_key if area_key is None else area_key
    eccentricity_key = ab.eccentricity_key if eccentricity_key is None else eccentricity_key
    col2adata(areas, data, area_key)
    col2adata(eccentricity, data, eccentricity_key)
    ab.stop_timer()  # write to obs, stop timer manually


@doc
def cell_co_occurrence(data: AnnData,
                       export_key: str = "cell_co_occurrence",
                       **kwargs):
    """The likelihood of two type of cells occur simultaneously in a ROI.

    Parameters
    ----------
    data : {adata}
    export_key : {export_key}
    **kwargs : {analysis_kwargs}

    """

    ab = AnalysisBase(data, display_name="Cell co-occurrence", export_key=export_key, **kwargs)
    ab.check_cell_type()
    df = ab.type_counter()
    df = df.T
    # normalize it using mean, greater than mean suggest it's occurrence
    df = ((df - df.mean()) / (df.max() - df.min()) > 0).astype(int)
    df = df.T
    # generate combination of cell types
    cell_comb = [i for i in combinations_with_replacement(df.columns, 2)]

    index = []
    values = []
    for c in cell_comb:
        c1 = c[0]
        c2 = c[1]
        # if two type of cells are all 1, the result is 1, if one is 0, the result is 0
        co_occur = (df[c1] * df[c2]).to_numpy()
        index.append((c1, c2))
        values.append(co_occur)
        if c1 != c2:
            index.append((c2, c1))
            values.append(co_occur)
    ab.result = pd.DataFrame(
        data=np.array(values).T,
        index=df.index,
        columns=pd.MultiIndex.from_tuples(index, names=['type1', 'type2']),
    )
