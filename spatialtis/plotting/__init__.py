from bokeh.io import output_notebook
from IPython.display import set_matplotlib_formats

from spatialtis.config import CONFIG

from ._bar_plot import stacked_bar
from ._cell_map import cell_map
from ._expression3d import expression_map
from ._heatmap_sns import heatmap
from ._stacked_kde_sns import stacked_kde
from ._violin_plot import violin_plot
from .palette import colorcycle, get_colors, get_linear_colors, view_palette
from .wrapper import (
    cc_interactions,
    cell_co_occurrence,
    cell_communities_graph,
    cell_components,
    cell_density,
    cell_morphology,
    cell_type_graph,
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
