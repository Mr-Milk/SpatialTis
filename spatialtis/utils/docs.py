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
    # keys
    key="The key stores the data for plotting in `AnnData.uns`",
    cell_type_key="The key to store cell types in `AnnData.obs` (Default: `spatialtis.CONFIG.CELL_TYPE_KEY`)",
    centroid_key="The key to store cell centroid in `AnnData.obs` (Default: `spatialtis.CONFIG.CENTROID_KEY`)",
    area_key="The key to store cell area in `AnnData.obs` (Default: `spatialtis.CONFIG.AREA_KEY`)",
    eccentricity_key="The key to store cell eccentricity in `AnnData.obs` "
    "(Default: `spatialtis.CONFIG.ECCENTRICITY_KEY`)",
    shape_key="The key to store cell shape in `AnnData.obs` (Default: `spatialtis.CONFIG.SHAPE_KEY`)",
    community_key="The key to store cell communities in `AnnData.obs`",
    neighbors_key="The key to store cell neighbors in `AnnData.obs`",
    marker_key="The key to store markers in `AnnData.var` (Default: `spatialtis.CONFIG.MARKER_KEY`)",
    layers_key="The layer in `AnnData` to perform analysis",
    # plot
    plot_options="Pass to :class:`spatialtis.plotting.abc.PlotBase`",
    adata_plotting="`AnnData` object for plotting",
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
    data_df="Data in `pandas.DataFrame` used for plotting",
)


def doc(obj):
    docstring = getdoc(obj)
    for param_name, content in PARAMETERS_DOCSTRINGS.items():
        pattern = re.compile(f"({{{param_name}}})")
        docstring = re.sub(pattern, content, docstring)
    obj.__doc__ = docstring
    return obj
