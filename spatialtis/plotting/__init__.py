from ..config import WORKING_ENV
from bokeh.io import output_notebook

from .palette import colorcycle, get_colors, get_linear_colors, view_palette
from ._bar_plot import stacked_bar
from ._cell_map import cell_map
from ._violin_plot import violin_plot
from .wrapper import cell_components, cell_density


if WORKING_ENV == 'jupyter':
    output_notebook(hide_banner=True)
elif WORKING_ENV == 'zeppelin':
    output_notebook(hide_banner=True, notebook_type='zepplin')
