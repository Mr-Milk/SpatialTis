from bokeh.io import output_notebook
from IPython.display import set_matplotlib_formats

from spatialtis.config import CONFIG

from ._bar_plot import stacked_bar
from ._cell_cell_interaction import cc_interactions
from ._cell_map import cell_map
from ._cell_map_echarts import cell_map_echarts
from ._community_graph import graph_plot, graph_plot_interactive
from ._dot_matrixplot import dot_matrix
from ._dotplot import dotplot
from ._expression_map import expression_map
from ._heatmap_sns import heatmap
from ._sankey import sankey
from ._triangle_dotplot import tri_dotplot
from ._violin_plot import violin_plot
from .palette import colorcycle, get_colors, get_linear_colors, view_palette
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

if CONFIG.WORKING_ENV == "jupyter":
    output_notebook(hide_banner=True)
elif CONFIG.WORKING_ENV == "zeppelin":
    output_notebook(hide_banner=True, notebook_type="zepplin")

# enable retina mode for all devices
set_matplotlib_formats("retina")
