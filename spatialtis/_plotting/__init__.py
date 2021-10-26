import matplotlib as mpl
from matplotlib import cycler

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

SPATIALTIS_STYLE = {
    "lines.linewidth": 2,
    "lines.markeredgecolor": "white",
    "lines.markeredgewidth": 1,
    "lines.markersize": 7,
    "patch.linewidth": 1,
    "patch.facecolor": "C0",
    "patch.edgecolor": "black",
    "text.color": "#0b0b0b",
    "axes.edgecolor": "#D0D0E0",
    "axes.grid": True,
    "axes.grid.axis": "both",
    "axes.titlesize": 10,
    "axes.labelsize": 10,
    "axes.labelcolor": "#0b0b0b",
    "axes.spines.left": False,
    "axes.spines.bottom": False,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.prop_cycle": cycler(
        "color",
        [
            "#00BEFF",
            "#D4CA3A",
            "#FF6DAE",
            "#67E1B5",
            "#EBACFA",
            "#9E9E9E",
            "#F1988E",
            "#5DB15A",
            "#E28544",
            "#52B8AA",
        ],
    ),
    "xtick.color": "#0b0b0b",
    "xtick.labelsize": 10,
    "ytick.color": "#0b0b0b",
    "ytick.labelsize": 10,
    "grid.color": "#93939c",
    "grid.linestyle": "--",
    "grid.alpha": 0.2,
    "figure.titlesize": 10,
    "figure.titleweight": "bold",
    "font.size": 10,
}

mpl.rcParams.update(SPATIALTIS_STYLE)