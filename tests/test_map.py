import matplotlib.pyplot as plt

import spatialtis as st
import spatialtis.plotting as sp
from spatialtis import Config

Config.reset()
Config.verbose = False


def test_cell_map(data2d, data2d_keys, data_wkt, data3d, data3d_keys):
    sp.cell_map(data2d, "ROI_1", roi_key="ROI")
    plt.close()
    sp.cell_map(data2d, "ROI_1", roi_key="ROI", cell_type_key="cell_type")
    plt.close()
    sp.cell_map(data2d_keys, "ROI_1", roi_key="ROI", cell_type_key="cell_type", centroid_key=('x', 'y'))
    plt.close()
    sp.cell_map(data_wkt, "ROI_1", roi_key="ROI", cell_type_key="cell_type")
    plt.close()
    sp.cell_map(data3d, "ROI_1", roi_key="ROI")
    plt.close()
    sp.cell_map(data3d, "ROI_1", roi_key="ROI", cell_type_key="cell_type")
    plt.close()
    sp.cell_map(data3d_keys, "ROI_1", roi_key="ROI", cell_type_key="cell_type", centroid_key=('x', 'y', 'z'))
    plt.close()


def test_exp_map(data2d, data2d_keys, data_wkt, data3d, data3d_keys):
    plt.close()
    sp.expression_map(data2d, "ROI_1", data2d.var.index[0], roi_key="ROI")
    plt.close()
    sp.expression_map(data2d, "ROI_1", data2d.var.index[0], roi_key="ROI", cell_type_key="cell_type")
    plt.close()

    sp.expression_map(data2d_keys, "ROI_1", data2d_keys.var.index[0],
                      roi_key="ROI", cell_type_key="cell_type", centroid_key=('x', 'y'))
    plt.close()

    sp.expression_map(data_wkt, "ROI_1", data_wkt.var.index[0], roi_key="ROI", cell_type_key="cell_type")
    plt.close()

    sp.expression_map(data3d, "ROI_1", data3d.var.index[0], roi_key="ROI")
    plt.close()
    sp.expression_map(data3d, "ROI_1", data3d.var.index[0], roi_key="ROI", cell_type_key="cell_type")
    plt.close()

    sp.expression_map(data3d_keys, "ROI_1", data3d_keys.var.index[0],
                      roi_key="ROI", cell_type_key="cell_type", centroid_key=('x', 'y', 'z'))
    plt.close()


def test_neighbors_map(data2d, data2d_keys, data_wkt):
    plt.close()
    st.find_neighbors(data2d, roi_key="ROI")
    sp.neighbors_map(data2d, "ROI_1", roi_key="ROI")
    plt.close()
    sp.neighbors_map(data2d, "ROI_1", roi_key="ROI", cell_type_key="cell_type")
    plt.close()

    st.find_neighbors(data2d_keys, roi_key="ROI", centroid_key=('x', 'y'))
    sp.neighbors_map(data2d_keys, "ROI_1", roi_key="ROI", cell_type_key="cell_type", centroid_key=('x', 'y'))
    plt.close()

    st.find_neighbors(data_wkt, roi_key="ROI")
    sp.neighbors_map(data_wkt, "ROI_1", roi_key="ROI", cell_type_key="cell_type")
    plt.close()
