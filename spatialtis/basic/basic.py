from itertools import combinations_with_replacement
from typing import Optional, List

import numpy as np
import pandas as pd
from anndata import AnnData
from spatialtis_core import multipoints_bbox, polygons_area, multipolygons_area

from spatialtis.abc import AnalysisBase
from spatialtis.utils import doc, read_points, read_shapes, col2adata_obs
from .utils import bbox_eccentricity


@doc
class cell_components(AnalysisBase):
    """Count the proportion of each types of cells in each group

    Args:
        data: {adata}
        **kwargs: {analysis_kwargs}

    """

    def __init__(self,
                 data: AnnData,
                 exp_obs: Optional[List[str]] = None,
                 roi_key: Optional[List[str]] = None,
                 cell_type_key: Optional[str] = None,
                 export_key: Optional[str] = None,
                 ):
        super().__init__(data,
                         exp_obs=exp_obs,
                         roi_key=roi_key,
                         cell_type_key=cell_type_key,
                         export_key=export_key)

        self.result = self.type_counter()


@doc
class cell_density(AnalysisBase):
    """Calculating cell density in each ROI

    The size of each ROI will be auto-computed, it's the area of convex hull of all the cells in a ROI

    Args:
        data: {adata}
        ratio: The ratio between the unit used in your dataset and real length unit, default is 1.0;
               ratio = Dataset unit / real length unit;
               For example, if the resolution of your dataset is 1μm, but you want to use 1mm as unit,
               then you should set the ratio as 0.001, 1 pixels represent 0.001mm length.
        **kwargs: {analysis_kwargs}

    """

    def __init__(self, data: AnnData, ratio: float = 1.0, **kwargs):
        super().__init__(data, **kwargs)
        df = self.type_counter()

        area = []
        for roi_name, roi_data in self.roi_iter():
            points = read_points(roi_data, self.centroid_key)
            area.append(polygons_area(points))

        area = np.asarray(area) * (ratio * ratio)
        self.results = df.div(area, axis=0)


@doc
class cell_morphology(AnalysisBase):
    """Cell morphology variation between different groups

    This function only works for data with cell shape information.
    The area is calculated using shoelace formula
    The eccentricity is assume that the cell is close to ellipse, the semi-minor and semi-major axis
    is get from the bbox side.

    Args:
        data: {adata}
        **kwargs: {analysis_kwargs}

    """

    def __init__(self, data: AnnData, metric: Optional[str] = None, **kwargs):
        super().__init__(data, **kwargs)

        shapes = read_shapes(self.data.obs, self.shape_key)
        areas = multipolygons_area(shapes)
        eccentricity = [bbox_eccentricity(bbox) for bbox in multipoints_bbox(shapes)]
        col2adata_obs(areas, self.data, self.area_key)
        col2adata_obs(eccentricity, self.data, self.eccentricity_key)
        self.stop_timer()  # write to obs, stop timer manually


@doc
class cell_co_occurrence(AnalysisBase):
    """The likelihood of two type of cells occur simultaneously in a ROI

    Args:
        data: {adata}
        **kwargs: {analysis_kwargs}

    """

    def __init__(self, data: AnnData, **kwargs):
        super().__init__(data, **kwargs)
        df = self.type_counter().T
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
            index.append((c1, c2))
            index.append((c2, c1))
            # if two type of cells are all 1, the result is 1, if one is 0, the result is 0
            co_occur = df[c1] * df[c2]
            values.append(co_occur)
            values.append(co_occur)

        self.result = pd.DataFrame(data=values, index=df.index, columns=pd.MultiIndex.from_tuples(index))
