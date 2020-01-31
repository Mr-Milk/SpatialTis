from ..config import WORKING_ENV
from bokeh.io import output_notebook

from .palette import colorcycle, get_colors, view_palette
from ._bar_plot import stacked_bar

if __name__ == "__main__":
    if WORKING_ENV == 'Jupyter':
        output_notebook(hide_banner=True)
    elif WORKING_ENV == 'Zeppelin':
        output_notebook(hide_banner=True, notebook_type='Zepplin')
