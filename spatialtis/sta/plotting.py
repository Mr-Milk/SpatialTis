from config import WORKING_ENV

from anndata import AnnData
import numpy as np
import pandas as pd

from typing import Optional


def cell_components(
        data: AnnData,
        count_base: str,
        metrics: str = 'percentage',
        read_key: str = 'cell_components',
        notebook: bool = True,
        save: Optional['png', 'svg', 'html'] = None
):
    if metrics == 'percentage':
        components = pd.DataFrame(index=all_results.index, columns=all_results.columns)
        for col, con in all_results.items():
            s = np.sum(con)
            components[col] = [count / s for count in con]
        container['data'] = components
