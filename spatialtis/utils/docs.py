import re
from inspect import getdoc

PARAMETERS_DOCSTRINGS = dict(
    adata="`AnnData` object to perform analysis",
    groupby="How to group your data in the plot",
    pval="The p-value threshold to determine significance",
    selected_types="Select your interested cell types",
    selected_markers="Select your interested markers",
    analysis_kwargs="Pass to :class:`spatialtis.abc.AnalysisBase`",
    tree_kwargs="The keyword arguments that pass to the boosting tree class, (Default: n_jobs=-1, random_state=0)",
    quad="A tuple (X, Y), Use a grid that is X * Y to tessellation your ROI",
    rect_side="A tuple (X, Y), Use the rectangle with X * Y side to tessellation your ROI",
    # keys
    export_key="The key used to store result",
    key="The key stores the data for plotting in `AnnData.uns`",
    exp_obs="A series of keys that store metadata in `AnnData.obs` (Default: :code:`spatialtis.Config.exp_obs`)",
    roi_key="The key that specify the ROI in `AnnData.obs` (Default: :code:`spatialtis.Config.roi_key`)",
    cell_type_key="The key to store cell types in `AnnData.obs` (Default: :code:`spatialtis.Config.cell_type_key`)",
    centroid_key="The key to store cell centroid in `AnnData.obs` (Default: :code:`spatialtis.Config.centroid_key`)",
    shape_key="The key to store cell shape in `AnnData.obs` (Default: :code:`spatialtis.Config.shape_key`)",
    community_key="The key to store cell communities in `AnnData.obs`",
    neighbors_key="The key to store cell neighbors in `AnnData.obs`",
    marker_key="The key to store markers in `AnnData.var` (Default: :code:`spatialtis.Config.marker_key`). "
               "If not specific, will use the index of `AnnData.var`",
    layer_key="The layer in `AnnData` to perform analysis",
    # plot
    agg="How to aggregate data, eg. sum, mean, median...",
    plot_options="Pass to :class:`spatialtis._plotting.abc.PlotBase`",
    adata_plotting="`AnnData` object for _plotting",
    roi="""A Dict use to select which ROI to display, 
             eg: :code:`{"Patients": "Patient 1", "ROI": "ROI3"}`, "Patients" and "ROI" are keys in `AnnData.obs`""",
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
