from typing import List, Optional

import pandas as pd
from anndata import AnnData

from spatialtis.config import ANALYSIS
from spatialtis.plotting.base import bar_static, violin_static
from spatialtis.utils import doc, get_result


@doc
def spatial_heterogeneity(
    data: AnnData,
    groupby: Optional[List[str]] = None,
    key: Optional[str] = None,
    metric: str = "heterogeneity",
    **kwargs,
):
    """Visualization for cell morphology

    Args:
        data: {adata_plotting}
        groupby: {groupby}
        key: {key}
        metric: "heterogeneity" or "KL", "KL" only available if you use shannon entropy
        **kwargs: Pass to :class:`spatialtis.plotting.base.violin_static` or
            :class:`spatialtis.plotting.base.bar_static`

    """
    if key is None:
        key = ANALYSIS["spatial_heterogeneity"].last_used_key

    if metric != "heterogeneity":
        metric = "KL"

    df, params = get_result(data, key, params=True)
    exp_obs = params["exp_obs"]

    if groupby is None:
        groupby = [exp_obs[0]]

    stacked = False
    names = []
    value = []
    for n, g in df.groupby(groupby):
        if isinstance(n, List):
            names.append("_".join(n))
        else:
            names.append(n)
        value.append(g[metric])
        if len(g) <= 1:
            stacked = True
        else:
            stacked = False

    if stacked:
        df = pd.DataFrame({"groups": names, metric: value})
        return bar_static(
            df, "groups", metric, **kwargs, saved_name="spatial_heterogeneity",
        )
    else:
        return violin_static(
            df,
            groupby,
            metric,
            hue=groupby[-1],
            **kwargs,
            saved_name="spatial_heterogeneity",
        )
