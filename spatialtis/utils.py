import matplotlib.pyplot as plt
import pandas as pd

from typing import Sequence


def df2adata_uns(df, adata, key, overwrite=False):
    """preserve all info in pd.DataFrame as dict, and write to anndata.uns

    """
    container = dict(
        df=str(df.to_dict()),
        iname=df.index.names,
        colname=list(df.columns.names),
    )

    keys = adata.uns.keys()
    if (key in keys) & (not overwrite):
        raise KeyError(f"{key} already exists, if you want to rewrite, set overwrite=True")

    adata.uns[key] = container

    print(
        f"""Finished!
    Add to AnnData object
    uns: '{key}' """
    )


def col2adata_obs(col: Sequence,
                  adata, key, overwrite=False):
    """preserve all info in pd.DataFrame as dict, and write to anndata.uns

    """
    keys = adata.uns.keys()
    if (key in keys) & (not overwrite):
        raise KeyError(f"{key} already exists, if you want to rewrite, set overwrite=True")

    adata.obs[key] = col

    print(
        f"""Finished!
    Add to AnnData object
    obs: '{key}' """
    )


def adata_uns2df(adata, key):

    container = adata.uns[key]
    df = pd.DataFrame(eval(container["df"]))
    df.index.set_names(container["iname"], inplace=True)
    df.columns.set_names(container["colname"], inplace=True)

    return df


def filter_adata(adata, groupby, type_col, *keys, selected_types=None, reset_index=True):

    keys = [k for k in keys if k is not None]
    df = adata.obs[groupby+keys+[type_col]]
    if selected_types is not None:
        df = df[df[type_col].isin(selected_types)]
    if reset_index:
        df = df.reset_index()

    return df


def plot_polygons(polygons):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    for i in polygons:
        ax.plot(
            *i.exterior.xy,
            color="#6699cc",
            alpha=0.7,
            linewidth=3,
            solid_capstyle="round",
            zorder=2,
        )
    ax.set_title("Polygons")
