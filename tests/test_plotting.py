import pandas as pd
import pytest

import spatialtis as st
import spatialtis._plotting as sp
from spatialtis import Config

plotting_kwargs = dict(
    title="just a title",
    xaxis_title="x title",
    yaxis_title="y title",
    legend_title="legend_title",
    xtickslabel_rotation=0,
    ytickslabel_rotation=0,
)

ROI = {"Patient": "HPAP005", "Part": "Tail", "ROI": "ROI1"}
TYPES = ["12", "14", "8"]
MARKERS = ['CD20', 'CD57', 'C-peptide', 'Ghrelin']
GORDER = {"Patient": ["HPAP005", "HPAP002"]}


def test_set_config(tmpdir):
    Config.exp_obs = ["Patient", "Part", "ROI"]
    Config.cell_type_key = "leiden"
    Config.MARKER_KEY = "Markers"
    Config.env = None
    Config.auto_save = True
    Config.SAVE_PATH = tmpdir


@pytest.mark.xfail
def test_colorcycle_failed():
    sp.colorcycle("V")


def test_get_colors():
    sp.get_colors(100, "Spectral")
    sp.get_colors(100, ["Spectral", "Set3"])


def test_get_linear_colors():
    sp.get_linear_colors("Spectral")
    sp.get_linear_colors(["Spectral", "Set3"])


def test_cell_map(data):
    sp.cell_map(data, ROI, use_shape=True, **plotting_kwargs)
    sp.cell_map(data, ROI, selected_types=TYPES)
    sp.cell_map(data, ROI, use_shape=True, use="interactive")
    sp.cell_map(data, ROI, selected_types=TYPES, use="interactive")


def test_expression_map(data):
    sp.expression_map(data, ROI, "CD20", **plotting_kwargs)
    sp.expression_map(data, ROI, "CD20", use="interactive").render()


def test_neighbors_map(data):
    sp.neighbors_map(data, ROI, **plotting_kwargs)
    sp.neighbors_map(data, ROI, use="interactive")


def test_community_map(data):
    sp.community_map(data, ROI, **plotting_kwargs)
    sp.community_map(data, ROI, use="interactive")


def test_cell_components(data):
    # one group
    sp.cell_components(data, ["Patient"], selected_types=TYPES, **plotting_kwargs)
    sp.cell_components(data, ["Patient"], direction="horizontal", group_order=GORDER)
    sp.cell_components(data, ["Patient", "Part"])
    sp.cell_components(data, ["Patient", "Part"], direction="horizontal", group_order=GORDER)

    sp.cell_components(data, ["Patient"], use="interactive", selected_types=TYPES, **plotting_kwargs)
    sp.cell_components(data, ["Patient"], use="interactive", direction="horizontal", group_order=GORDER)
    sp.cell_components(data, ["Patient", "Part"], use="interactive", selected_types=TYPES)
    sp.cell_components(data, ["Patient", "Part"], use="interactive", direction="horizontal", group_order=GORDER)


def test_cell_density(data):
    sp.cell_density(data, ["Patient"], selected_types=TYPES, **plotting_kwargs)
    sp.cell_density(data, ["Patient"], direction="horizontal", group_order=GORDER)
    sp.cell_density(data, ["Patient", "Part"])
    sp.cell_density(data, ["Patient", "Part"], direction="horizontal", group_order=GORDER)


def test_cell_morphology(data):
    sp.cell_morphology(data, ['Patient', 'Part'], selected_types=TYPES)


def test_cell_co_occurrence(data):
    sp.cell_co_occurrence(data, selected_types=TYPES)
    sp.cell_co_occurrence(data, use='heatmap')


def test_spatial_distribution(data):
    # sp.spatial_distribution(data, selected_types=TYPES)
    sp.spatial_distribution(data, use='heatmap')


def test_spatial_heterogeneity(data):
    sp.spatial_heterogeneity(data)


def test_neighborhood_analysis(data):
    sp.neighborhood_analysis(data, selected_types=TYPES, use="heatmap")
    sp.neighborhood_analysis(data, use="dot_matrix")
    sp.neighborhood_analysis(data, use="graph_static")
    sp.neighborhood_analysis(data, use="graph_interactive").render()


def test_neighborhood_analysis_order(data):
    sp.neighborhood_analysis(data, selected_types=TYPES, use="heatmap", key="neighborhood_order")
    sp.neighborhood_analysis(data, use="dot_matrix", key="neighborhood_order")
    sp.neighborhood_analysis(data, use="graph_static", key="neighborhood_order")
    sp.neighborhood_analysis(data, use="graph_interactive", key="neighborhood_order").render()


def test_spatial_enrichment(data):
    sp.spatial_enrichment_analysis(data, selected_markers=MARKERS, use="heatmap")
    sp.spatial_enrichment_analysis(data, use="dot")
    sp.spatial_enrichment_analysis(data, use="graph_static")
    sp.spatial_enrichment_analysis(data, use="graph_interactive").render()


def test_spatial_enrichment_order(data):
    sp.spatial_enrichment_analysis(data, selected_markers=MARKERS, use="heatmap", key="enrichment_order")
    sp.spatial_enrichment_analysis(data, use="dot", key="enrichment_order")
    sp.spatial_enrichment_analysis(data, use="graph_static", key="enrichment_order")
    sp.spatial_enrichment_analysis(data, use="graph_interactive", key="enrichment_order").render()


def test_co_expression(data):
    test_df = pd.DataFrame(data=[['CD20', 'CD20', 0.5, 0.0001],
                                 ['CD57', 'CD57', 0.5, 0.0001],
                                 ['CD20', 'CD57', 0.5, 0.0001]],
                           columns=['marker1', 'marker2', 'corr', 'pvalue'])
    st.utils.df2adata_uns(test_df, data, 'co_expression')
    sp.spatial_co_expression(data, selected_markers=MARKERS, use="heatmap")
    sp.spatial_co_expression(data, use="graph_static")
    sp.spatial_co_expression(data, use="graph_interactive")


def test_NCDMarkers(data):
    test_df = pd.DataFrame(data=[['CD20', 'B', 0.5, 0.4, 0.0001],
                                 ['CD23', 'C', 0.9, -0.4, 0.0001]],
                           columns=['marker', 'neighbor_type', 'dependency', 'corr', 'pvalue'])
    st.utils.df2adata_uns(test_df, data, 'ncd_markers')
    sp.NCDMarkers(data)


def test_NMDMarkers(data):
    test_df = pd.DataFrame(data=[['CD20', 'B', 0.5, 0.4],
                                 ['CD23', 'C', 0.9, -0.4]],
                           columns=["neighbor_marker", "marker", "dependency", "corr"])
    st.utils.df2adata_uns(test_df, data, 'nmd_markers')
    sp.NMDMarkers(data)


def test_base_bar():
    df = pd.DataFrame({"x": list(range(10)),
                       "y": list(range(10)),
                       })
    sp.base.bar_static(df, x="x", y="y")
    sp.base.bar_static(df, x="x", y="y", direction="horizontal")
