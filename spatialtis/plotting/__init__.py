from bokeh.io import output_notebook
# enable retina mode for all devices
from IPython.display import set_matplotlib_formats

from spatialtis.config import CONFIG

from ._bar_plot import stacked_bar
from ._cell_map import cell_map
from ._heatmap_sns import heatmap
from ._stacked_kde_sns import stacked_kde
from ._violin_plot import violin_plot
from .palette import colorcycle, get_colors, get_linear_colors, view_palette
from .wrapper import (cell_co_occurrence, cell_components, cell_density,
                      cell_morphology, neighborhood_analysis,
                      spatial_distribution, spatial_enrichment_analysis,
                      spatial_heterogeneity)

if CONFIG.WORKING_ENV == "jupyter":
    output_notebook(hide_banner=True)
elif CONFIG.WORKING_ENV == "zeppelin":
    output_notebook(hide_banner=True, notebook_type="zepplin")


set_matplotlib_formats("retina")
