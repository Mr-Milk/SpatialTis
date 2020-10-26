from .base import (
    DotMatrix,
    TriDotMatrix,
    colorcycle,
    dotplot,
    get_colors,
    get_linear_colors,
    graph_plot,
    graph_plot_interactive,
    heatmap,
    sankey,
    stacked_bar,
    tri_dotplot,
    violin_plot,
)
from .roi_viz import cell_communities, cell_map, cell_neighbors, expression_map
from .wrapper import (
    cell_co_occurrence,
    cell_components,
    cell_density,
    cell_morphology,
    exp_neighcells,
    neighborhood_analysis,
    spatial_distribution,
    spatial_enrichment_analysis,
    spatial_heterogeneity,
)

# enable retina mode for all devices
try:
    from IPython.display import set_matplotlib_formats

    set_matplotlib_formats("retina")
except ImportError:
    pass
