import spatialtis as st
import spatialtis.plotting as sp
from spatialtis import Config

type_order = ["a", "b", "c"]


# Basic analysis
def test_2d(data2d, data2d_keys, data_wkt):
    Config.loads(data2d)
    st.cell_components(data2d)
    st.cell_density(data2d)
    st.cell_co_occurrence(data2d)

    st.cell_density(data2d_keys, centroid_key=["x", "y"])
    st.cell_density(data_wkt)


# Not supported yet
# def test_3d(data3d):
#     st.cell_density(data3d)


def test_basic_plot(data2d):
    Config.loads(data2d)
    sp.cell_components(data2d)
    sp.cell_components(data2d, orient="h", type_order=type_order)

    sp.cell_density(data2d)
    sp.cell_density(data2d, groupby="ROI", type_order=type_order)

    sp.cell_co_occurrence(data2d)
    sp.cell_co_occurrence(data2d, groupby=["ROI"], type_order=type_order)


def test_morphology(data_shape):
    Config.loads(data_shape)
    st.cell_morphology(data_shape)
