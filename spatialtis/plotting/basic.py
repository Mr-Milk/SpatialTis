from typing import Optional, List

import pandas as pd
import seaborn as sns
from anndata import AnnData
from milkviz import stacked_bar, dot, anno_clustermap
from milkviz.utils import mask_triu

from spatialtis.utils import doc, get_result
from .utils import pairs_to_adj
from .. import Config


@doc
def cell_components(data: AnnData,
                    groupby: Optional[str] = None,
                    key: str = "cell_components",
                    orient: str = "v",
                    type_order: Optional[List[str]] = None,
                    **plot_options,
                    ):
    """Visualize cell components result

    Args:
        data:
        groupby:
        key:
        orient:
        type_order:
        **plot_options:

    Returns:

    """

    data = get_result(data, key)
    groupby = data.index.names[0] if groupby is None else groupby
    if type_order is not None:
        data = data.loc[:, type_order]
    data = data.groupby(groupby).sum().melt(ignore_index=False, value_name="Count").reset_index()
    data = data.rename(columns={'cell type': 'Cell Type'})
    if orient == "v":
        return stacked_bar(data,
                           x=groupby,
                           y="Count",
                           stacked="Cell Type",
                           orient=orient,
                           **plot_options)
    else:
        return stacked_bar(data,
                           y=groupby,
                           x="Count",
                           stacked="Cell Type",
                           orient=orient,
                           **plot_options)


@doc
def cell_density(data: AnnData,
                 groupby: Optional[str] = None,
                 key: str = "cell_density",
                 type_order: Optional[List[str]] = None,
                 **plot_options,
                 ):
    data = get_result(data, key)
    if type_order is not None:
        data = data.loc[:, type_order]
    data = pd.melt(data, ignore_index=False, value_name="density").reset_index()
    options = dict(palette="tab20", **plot_options)
    if groupby is None:
        ax = sns.boxplot(data=data, x="cell type", y="density", **options)
    else:
        ax = sns.boxplot(data=data, x=groupby, y="density", hue="cell type", **options)
        ax.legend(loc="upper left", bbox_to_anchor=(1.05, 0, 1, 1), title="cell type", frameon=False)
    ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
    ax.set(xlabel="Cell Type", ylabel="Density")
    return ax


@doc
def cell_morphology(data: AnnData,
                    groupby: Optional[str] = None,
                    key: str = "area",
                    type_order: Optional[List[str]] = None,
                    cell_type_key: Optional[str] = None,
                    **plot_options,
                    ):
    cell_type_key = Config.cell_type_key if cell_type_key is None else cell_type_key
    options = dict(palette="tab20", **plot_options)
    query_keys = [key, cell_type_key] if groupby is None else [key, cell_type_key, groupby]
    pdata = data.obs[query_keys]
    if groupby is None:
        ax = sns.boxplot(data=pdata, x=cell_type_key, y=key, order=type_order, **options)
    else:
        ax = sns.boxplot(data=pdata, x=groupby, y=key, hue=cell_type_key, hue_order=type_order, **options)
        ax.legend(loc="upper left", bbox_to_anchor=(1.05, 0, 1, 1), title="cell type", frameon=False)
    ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
    ax.set(xlabel="Cell Type", ylabel=key.capitalize())
    return ax


@doc
def cell_co_occurrence(data: AnnData,
                       use: str = "dot",  # dot, heatmap
                       groupby: Optional[List] = None,
                       key: str = "cell_co_occurrence",
                       type_order: Optional[List[str]] = None,
                       order: bool = True,
                       **plot_options,
                       ):
    data = get_result(data, key)
    if use == "dot":
        data = pd.DataFrame(data.sum().reset_index())
        data = pairs_to_adj(data, type_order)
        plot_data = data.to_numpy()
        xticklabels = data.columns
        yticklabels = data.index.tolist()
        if not order:
            plot_data = mask_triu(plot_data)
            xticklabels = xticklabels[::-1]
        plot_options = {"dot_patch": "pie", **plot_options}
        return dot(plot_data,
                   xticklabels=xticklabels,
                   yticklabels=yticklabels,
                   **plot_options
                   )
    else:
        groupby = data.index.names[0] if groupby is None else groupby
        plot_options = {"categorical_cbar": ['Non-Occur', 'Co-occur'],
                        "col_legend_split": False,
                        "col_legend_title": "Cell Type",
                        "cbar_title": "Occurrence",
                        **plot_options}
        return anno_clustermap(data,
                               col_colors=['type1', 'type2'],
                               row_colors=groupby,
                               **plot_options
                               )
