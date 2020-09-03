from ._bar_plot import stacked_bar
from ._cell_map import cell_map
from ._community_graph import graph_plot, graph_plot_interactive
from ._dot_matrix import DotMatrix
from ._dotplot import dotplot
from ._expression_map import expression_map
from ._heatmap_sns import heatmap
from ._sankey import sankey
from ._tri_dot_matrix import TriDotMatrix
from ._triangle_dotplot import tri_dotplot
from ._violin_plot import violin_plot
from .palette import colorcycle, get_colors, get_linear_colors
from .wrapper import (
    cell_co_occurrence,
    cell_communities,
    cell_components,
    cell_density,
    cell_morphology,
    cell_neighbors,
    exp_neighcells,
    neighborhood_analysis,
    spatial_distribution,
    spatial_enrichment_analysis,
    spatial_heterogeneity,
)

"""
# enable retina mode for all devices
try:
    from IPython.display import set_matplotlib_formats

    set_matplotlib_formats("retina")
except ImportError:
    pass
"""
