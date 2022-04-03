import spatialtis as st
import spatialtis.plotting as sp
from spatialtis import Config

Config.reset()

Config.verbose = False
Config.roi_key = "ROI"
Config.cell_type_key = "cell_type"
type_order = ["a", "b", "c"]


def test_cell_dispersion(data2d):
    st.cell_dispersion(data2d, r=50, method="id")
    st.cell_dispersion(data2d, quad=(10, 10), method="morisita")
    st.cell_dispersion(data2d, method="clark_evans")


def test_cell_dispersion_3d(data3d):
    st.cell_dispersion(data3d, r=50, method="id")


def test_cell_dispersion_plot(data2d):
    sp.cell_dispersion(data2d, type_order=type_order)
    sp.cell_dispersion(data2d, use="heatmap", type_order=type_order)


def test_spatial_hetero(data2d):
    st.spatial_heterogeneity(data2d, method="shannon")
    st.spatial_heterogeneity(data2d, method="leibovici", d=10)
    st.spatial_heterogeneity(data2d, method="altieri", cut=3)


def test_spatial_hetero_3d(data3d):
    st.spatial_heterogeneity(data3d, method="shannon")
    st.spatial_heterogeneity(data3d, method="leibovici", d=10)
    st.spatial_heterogeneity(data3d, method="altieri", cut=3)


def test_spatial_hetero_plot(data2d):
    sp.spatial_heterogeneity(data2d)


def test_hotspot(data2d):
    st.hotspot(data2d, rect_side=(10, 10))

