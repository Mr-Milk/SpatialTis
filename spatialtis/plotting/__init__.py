# from .base import (
#     DotMatrix,
#     TriDotMatrix,
#     colorcycle,
#     dotplot,
#     get_colors,
#     get_linear_colors,
#     graph_plot,
#     graph_plot_interactive,
#     heatmap,
#     sankey,
#     stacked_bar,
#     tri_dotplot,
#     violin_plot,
# )
import matplotlib as mpl

from .api import (
    NCDMarkers,
    NMDMarkers,
    cell_co_occurrence,
    cell_components,
    cell_density,
    cell_map,
    cell_morphology,
    community_map,
    expression_map,
    neighborhood_analysis,
    neighbors_map,
    spatial_co_expression,
    spatial_distribution,
    spatial_enrichment_analysis,
    spatial_heterogeneity,
)
from .base import get_colors, get_linear_colors

try:
    import themepy

    theme = themepy.Theme()
    theme.set_theme("gadfly")
except ImportError:
    pass

mpl.rcParams["font.size"] = 8

# enable retina mode for all devices
try:
    from IPython.display import set_matplotlib_formats

    set_matplotlib_formats("retina")
except ImportError:
    pass
