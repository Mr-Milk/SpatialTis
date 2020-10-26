import re
from inspect import getdoc

PARAMETERS_DOCSTRINGS = dict(
    adata="`AnnData` object to perform analysis",
    adata_plotting="`AnnData` object for plotting",
    groupby="How your experiments data grouped, (Default: `spatialtis.CONFIG.EXP_OBS`)",
    pval="The p-value threshold to determine significance",
    selected_types="Select interested cell types to perform analysis",
    selected_markers="Select interested markers to perform analysis",
    n="A `spatialtis.Neighbors` object, neighbors are already computed",
    type_key="The key to store cell types in `AnnData.obs` (Default: `spatialtis.CONFIG.CELL_TYPE_KEY`)",
    centroid_key="The key to store cell centroid in `AnnData.obs` (Default: `spatialtis.CONFIG.CENTROID_KEY`)",
    shape_key="The key to store cell shape in `AnnData.obs` (Default: `spatialtis.CONFIG.SHAPE_KEY`)",
    community_key="The key to store cell communities in `AnnData.obs`",
    neighbors_key="The key to store cell neighbors in `AnnData.obs`",
    marker_key="The key to store markers in `AnnData.var` (Default: `spatialtis.CONFIG.MARKER_KEY`)",
    layers_key="The layer in anndata to perform analysis",
    export="Whether export the result to `AnnData.uns`",
    export_key="The name of key used to stored the exported result",
    return_df="Whether to return the result dataframe",
    mp="Whether to enable parallel processing (Default: `spatialtis.CONFIG.MULTI_PROCESSING`)",
    key="The key stores the data for plotting in `AnnData.uns`",
    query="""A Dict use to select which ROI to display, 
             eg: {"Patients": "Patient 1", "ROI": "ROI3"}, "Patients" and "ROI" are keys in `AnnData.obs`""",
    group_order="""A Dict to control the order of axis-label
                    eg: {'stage':['stage 1', 'stage 2'],'case':[1,3,6]}
                    The plot will display the level of stage following ['stage 1', 'stage 2'] and
                    the level of case following [1,3,6]""",
    direction="Options are 'vertical' and 'horizontal'",
    xaxis_title="The title of X axis",
    yaxis_title="The title of Y axis",
    xlabels="The labels marked on X axis",
    ylabels="The labels marked on Y axis",
    xlabel_rotation="Number of degree to rotate X-axis's label",
    ylabel_rotation="Number of degree to rotate Y-axis's label",
    legend_title="The title of legend",
    annotate="Whether to show value number",
    alpha="Alpha value for opacity",
    renderer="""Options are "canvas" and "svg"; "svg" is not perfect at this moment.""",
    theme="Please go to `theme <https://pyecharts.org/#/zh-cn/themes>`_ for details",
    title="The title of the plot",
    size="The size of the plot in pixels",
    palette="""Control the color of plot, sequence of color in hex, or
            `name of palettes <https://docs.bokeh.org/en/latest/docs/reference/palettes.html>`_""",
    display="Whether to display the plot",
    save="The path to save your plot",
    return_plot="Whether to return the plot instance",
    show_ticks="Whether to show the minor ticks",
)


def reuse_docstring():
    def wrapper(func):
        docstring = getdoc(func)
        for param_name, content in PARAMETERS_DOCSTRINGS.items():
            pattern = re.compile(f"({{{param_name}}})")
            docstring = re.sub(pattern, content, docstring)
        func.__doc__ = docstring
        return func

    return wrapper
