from ast import literal_eval
from itertools import combinations_with_replacement
from typing import Optional

import numpy as np
import pandas as pd
from anndata import AnnData
from shapely.geometry import MultiPoint

from spatialtis.abc import AnalysisBase
from spatialtis.typing import Number
from spatialtis.utils import doc


@doc
class cell_components(AnalysisBase):
    """Count the proportion of each types of cells in each group

    Args:
        data: {adata}
        **kwargs: {analysis_kwargs}

    """

    def __init__(self, data: AnnData, **kwargs):
        super().__init__(data, task_name="cell_components", **kwargs)
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

    def __init__(self, data: AnnData, ratio: Number = 1.0, **kwargs):
        super().__init__(data, task_name="cell_density", **kwargs)
        counter = self.type_counter()
        df = counter.iloc[:, len(self.exp_obs) : :]

        groups = data.obs[self.exp_obs + [self.centroid_key]].groupby(self.exp_obs)
        area = [
            MultiPoint(
                [literal_eval(c) for c in g[self.centroid_key].tolist()]
            ).convex_hull.area
            for _, g in groups
        ]
        area = np.asarray(area) * ratio * ratio
        results = df.div(area, axis=0)
        results = pd.concat([counter.loc[:, self.exp_obs], results], axis=1)
        self.result = pd.melt(results, id_vars=self.exp_obs, var_name="type")


@doc
class cell_morphology(AnalysisBase):
    """Cell morphology variation between different groups

    Args:
        data: {adata}
        metric_key: Which key in `AnnData.obs` to measure morphology (Default: eccentricity key)
        **kwargs: {analysis_kwargs}

    """

    def __init__(self, data: AnnData, metric_key: Optional[str] = None, **kwargs):
        super().__init__(data, task_name="cell_morphology", **kwargs)
        if metric_key is None:
            metric_key = self.eccentricity_key
        result = data.obs.loc[:, self.exp_obs + [self.cell_type_key, metric_key]]
        result.rename(
            columns={self.cell_type_key: "type", metric_key: "value"}, inplace=True
        )
        self.result = result


@doc
class cell_co_occurrence(AnalysisBase):
    """The likelihood of two type of cells occur simultaneously in a ROI

    Args:
        data: {adata}
        **kwargs: {analysis_kwargs}

    """

    def __init__(self, data: AnnData, **kwargs):
        super().__init__(data, task_name="cell_co_occurrence", **kwargs)
        counter = self.type_counter()
        df = counter.iloc[:, len(self.exp_obs) : :]
        # normalize it using mean, greater than mean suggest it's occurrence
        df = ((df - df.mean()) / (df.max() - df.min()) > 0).astype(int)
        # generate combination of cell types
        cell_comb = [i for i in combinations_with_replacement(df.columns, 2)]
        chunks = []
        for c in cell_comb:
            # if two type of cells are all 1, the result is 1, if one is 0, the result is 0
            co = pd.DataFrame({"co_occur": df[c[0]] * df[c[1]]})
            co["type1"] = c[0]
            co["type2"] = c[1]
            co = pd.concat([co, counter.loc[:, self.exp_obs]], axis=1)
            chunks.append(co)
        self.result = pd.concat(chunks)[
            self.exp_obs + ["type1", "type2", "co_occur"]
        ].reset_index(drop=True)
