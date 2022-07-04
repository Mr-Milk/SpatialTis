import pandas as pd
from typing import List, Optional

COLOR_POOL = [
    "#1f77b4ff",
    "#aec7e8ff",
    "#ff7f0eff",
    "#ffbb78ff",
    "#2ca02cff",
    "#98df8aff",
    "#d62728ff",
    "#ff9896ff",
    "#9467bdff",
    "#c5b0d5ff",
    "#8c564bff",
    "#c49c94ff",
    "#e377c2ff",
    "#f7b6d2ff",
    "#bcbd22ff",
    "#dbdb8dff",
    "#17becfff",
    "#9edae5ff",
    "#e41a1cff",
    "#377eb8ff",
    "#4daf4aff",
    "#984ea3ff",
    "#ff7f00ff",
    "#ffff33ff",
    "#a65628ff",
    "#f781bfff",
]


def pairs_to_adj(data: pd.DataFrame, type_order: Optional[List] = None):
    data.columns = ["type1", "type2", "value"]
    data = pd.pivot(data, columns="type1", index="type2", values="value")
    if type_order is not None:
        data = data.loc[type_order, type_order]
    return data
