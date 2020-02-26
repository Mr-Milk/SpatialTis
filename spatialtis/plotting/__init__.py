from bokeh.io import output_notebook

from ..config import WORKING_ENV
from ._bar_plot import stacked_bar
from ._cell_map import cell_map
from ._violin_plot import violin_plot
from ._heatmap_sns import heatmap
from ._stacked_kde_sns import stacked_kde
from .palette import colorcycle, get_colors, get_linear_colors, view_palette
from .wrapper import cell_components, cell_density, cell_co_occurrence

if WORKING_ENV == "jupyter":
    output_notebook(hide_banner=True)
elif WORKING_ENV == "zeppelin":
    output_notebook(hide_banner=True, notebook_type="zepplin")

# enable retina mode for all devices
from IPython.display import set_matplotlib_formats
set_matplotlib_formats('retina')
