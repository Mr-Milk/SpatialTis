import spatialtis as st
import spatialtis.plotting as sp
from spatialtis import Config

Config.reset()

Config.verbose = False
Config.roi_key = "ROI"
Config.cell_type_key = "cell_type"
type_order = ["a", "b", "c"]


def test_find_neighbors(data2d, data3d, data_shape):
    st.find_neighbors(data2d, method="delaunay")
    st.find_neighbors(data2d, r=30, k=5, method="kdtree")

    st.find_neighbors(data3d, r=30, k=5, method="kdtree")

    st.find_neighbors(data_shape, scale=1.2, method="rtree", shape_key="cell_shape")


def test_spatial_autocorr(data2d, data3d):
    st.spatial_autocorr(data2d, method="moran_i")
    st.spatial_autocorr(data2d, method="geary_c")

    st.spatial_autocorr(data3d, method="moran_i")
    st.spatial_autocorr(data3d, method="geary_c")


def test_neighborhood(data2d, data3d):
    st.cell_interaction(data2d)
    st.cell_interaction(data3d, method="zscore")


def test_neighborhood_plot(data2d):
    sp.cell_interaction(data2d)
    sp.cell_interaction(data2d, use="heatmap", type_order=type_order)
    sp.cell_interaction(data2d, use="heatmap", plot_value="value")


def test_spatial_enrichment(data2d, data3d):
    st.spatial_enrichment(data2d)
    st.spatial_enrichment(data3d, threshold=0.1)


def test_spatial_enrichment_plot(data2d):
    sp.spatial_enrichment(data2d)


def test_spatial_coexp(data2d, data3d):
    st.spatial_coexp(data2d)
    st.spatial_coexp(data2d, use_cell_type=True)

    st.spatial_coexp(data3d)
    st.spatial_coexp(data3d, use_cell_type=True)


def test_ncd_markers(data2d, data3d):
    st.NCD_marker(data2d)
    st.NCD_marker(data3d)


def test_nmd_markers(data2d, data3d):
    st.NMD_marker(data2d)
    st.NMD_marker(data3d)
