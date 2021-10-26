from itertools import combinations_with_replacement
from typing import List, Optional

import numpy as np
import pandas as pd
from anndata import AnnData
from spatialtis_core import multipoints_bbox, multipolygons_area, polygons_area

from spatialtis.abc import AnalysisBase
from spatialtis.utils import col2adata_obs, doc, read_points, read_shapes
from .utils import bbox_eccentricity


@doc
class cell_components(AnalysisBase):
    """Count the proportion of each types of cells in each group

    Args:
        data: {adata}
        export_key: {export_key}
        **kwargs: {analysis_kwargs}

    """

    def __init__(
        self,
        data: AnnData,
        export_key: Optional[str] = None,
        **kwargs,
    ):
        super().__init__(data, export_key=export_key, **kwargs)

        result = self.type_counter()
        result.columns.name = 'cell type'
        self.result = result


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
        export_key: {export_key}
        **kwargs: {analysis_kwargs}

    """

    def __init__(self,
                 data: AnnData,
                 ratio: float = 1.0,
                 export_key: Optional[str] = None,
                 **kwargs
                 ):
        super().__init__(data, export_key=export_key, **kwargs)
        df = self.type_counter()

        area = []
        for roi_name, roi_data in self.roi_iter():
            points = read_points(roi_data, self.centroid_key)
            area.append(polygons_area(points))

        area = np.asarray(area) * (ratio * ratio)
        result = df.div(area, axis=0)
        result.columns.name = 'cell type'
        self.result = result


@doc
class cell_morphology(AnalysisBase):
    """Cell morphology variation between different groups

    This function only works for data with cell shape information.
    The area is calculated using shoelace formula
    The eccentricity is assume that the cell is close to ellipse, the semi-minor and semi-major axis
    is get from the bbox side.

    Args:
        data: {adata}
        area_key: The key to store cell area, Default: 'area'
        eccentricity_key: The key to store cell eccentricity, Default: 'eccentricity'
        **kwargs: {analysis_kwargs}

    """

    def __init__(self,
                 data: AnnData,
                 area_key: Optional[str] = None,
                 eccentricity_key: Optional[str] = None,
                 **kwargs):
        super().__init__(data, **kwargs)

        shapes = read_shapes(self.data.obs, self.shape_key)
        areas = multipolygons_area(shapes)
        eccentricity = [bbox_eccentricity(bbox) for bbox in multipoints_bbox(shapes)]
        area_key = self.area_key if area_key is None else area_key
        eccentricity_key = self.eccentricity_key if eccentricity_key is None else eccentricity_key
        col2adata_obs(areas, self.data, area_key)
        col2adata_obs(eccentricity, self.data, eccentricity_key)
        self.stop_timer()  # write to obs, stop timer manually


@doc
class cell_co_occurrence(AnalysisBase):
    """The likelihood of two type of cells occur simultaneously in a ROI

    Args:
        data: {adata}
        export_key: {export_key}
        **kwargs: {analysis_kwargs}

    """

    def __init__(self,
                 data: AnnData,
                 export_key: Optional[str] = None,
                 **kwargs):
        super().__init__(data, export_key=export_key, **kwargs)
        df = self.type_counter()
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
        self.result = pd.DataFrame(
            data=np.array(values).T,
            index=df.index,
            columns=pd.MultiIndex.from_tuples(index, names=['type1', 'type2']),
        )
