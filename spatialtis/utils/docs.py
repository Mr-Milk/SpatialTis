import re
from inspect import getdoc


def docstring(types, content):
    return f"""{types}
        {content}
    """


PARAMETERS_DOCSTRINGS = dict(
    adata=docstring("AnnData", "The `AnnData` to work with."),
    export_key=docstring("str", "The key used to store result."),
    analysis_kwargs="Config for the analysis, for details check :class:`spatialtis.abc.AnalysisBase`.",
    pval=docstring("float", "The p-value threshold to determine significance"),
    selected_types=docstring("list of str", "Select your interested cell types"),
    selected_markers=docstring("list of str", "Select your interested markers"),

    # keys
    exp_obs=docstring("str or list of str, default: :code:`Config.exp_obs`",
                      "The columns in `.obs` that tells how your experiments organized. "
                      "This will temporarily overwrite the global config."),
    roi_key=docstring("str, default: :code:`Config.roi_key`",
                      "The columns in `.obs` that specific ROI."
                      "This will temporarily overwrite the global config."),
    cell_type_key=docstring("str, default: :code:`Config.cell_type_key`",
                            "The column in `.obs` to store cell types."
                            "This will temporarily overwrite the global config."),
    centroid_key=docstring("str or list of str, default: :code:`Config.centroid_key`",
                           "The column in `.obs` or `.obsm` to store centroid. "
                           "This will temporarily overwrite the global config."),
    marker_key=docstring("str, default: :code:`Config.marker_key`",
                         "The column in `.var` to store marker names."
                         "By default will read from the index of `.var`. "
                         "This will temporarily overwrite the global config."),
    shape_key=docstring("str, default: :code:`Config.shape_key`",
                        "The column in `.var` to store cell shape in wkt format. "
                        "This will temporarily overwrite the global config."),
    layer_key=docstring("str", "The layer in `AnnData` to perform analysis"),
    # plotting
    adata_plotting=docstring("AnnData", "The `AnnData` object for plotting"),
    plot_key=docstring("str", "The key that store the analysis results"),
    groupby=docstring("str", "The key in :code:`Config.exp_obs` to stratify the data."),

    tree_kwargs="The keyword arguments that pass to the boosting tree class, (Default: n_jobs=-1, random_state=0)",

    community_key="The key to store cell communities in `AnnData.obs`",
    neighbors_key="The key to store cell neighbors in `AnnData.obs`",

    # plot
    agg="How to aggregate data, eg. sum, mean, median...",

    pyecharts_tips=""".. note::
    If using interactive plot here, you need to call `.render()` to display the plot in the notebook""",
    group_order="""A Dict to control the order of axis-label
                    eg: :code:`{'stage':['stage 1', 'stage 2'],'case':[1,3,6]}`
                    The plot will display the level of stage following :code:`['stage 1', 'stage 2']` and
                    the level of case following :code:`[1,3,6]`""",
    direction='"vertical" and "horizontal" (Default: "vertical")',
    annotate="Whether to show value number",
    alpha="Alpha value for opacity",
    data_df="Data in `pandas.DataFrame` used for _plotting",
)


def doc(obj):
    docstring = getdoc(obj)
    if docstring is not None:
        for param_name, content in PARAMETERS_DOCSTRINGS.items():
            pattern = re.compile(f"({{{param_name}}})")
            docstring = re.sub(pattern, content, docstring)
        obj.__doc__ = docstring
    return obj
