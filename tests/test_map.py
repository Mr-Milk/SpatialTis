import matplotlib.pyplot as plt

import spatialtis as st
import spatialtis.plotting as sp
from spatialtis import Config

Config.reset()
Config.verbose = False

ROIs = ["ROI_1", "ROI_2"]
Shape_ROIs = ["1", "2"]


def test_cell_map2d(data2d):
    Config.loads(data2d)
    sp.cell_map(data2d, ROIs)
    st.find_neighbors(data2d)
    sp.cell_map(data2d, ROIs, show_neighbors=True)
    sp.cell_map(data2d, ROIs, selected_types=['a', 'b', 'c'], show_neighbors=True)


def test_cell_map_shape(data_shape):
    Config.loads(data_shape)
    sp.cell_map(data_shape, Shape_ROIs, use_shape=True)
    sp.cell_map(data_shape, Shape_ROIs, use_shape=True, selected_types=["1", "2", "3"])


def test_cell_map3d(data3d, data3d_keys):
    Config.loads(data3d)
    sp.cell_map(data3d, ROIs)
    sp.cell_map(data3d, ROIs, selected_types=['a', 'b', 'c'])


def test_exp_map2d(data2d):
    Config.loads(data2d)
    markers = data2d.var.index[[0, 1, 2, 3]]
    sp.expression_map(data2d, ROIs, markers)
    sp.expression_map(data2d, ROIs, markers, x_axis="roi")
    sp.expression_map(data2d, ROIs, markers, selected_types=['a', 'b', 'c'])


def test_exp_map_shape(data_shape):
    Config.loads(data_shape)
    markers = data_shape.var['Markers'][[0, 1, 2, 3]]
    sp.expression_map(data_shape, Shape_ROIs, markers)
    sp.expression_map(data_shape, Shape_ROIs, markers, selected_types=["1", "2", "3"])


def test_exp_map3d(data3d):
    Config.loads(data3d)
    markers = data3d.var.index[[0, 1, 2, 3]]
    sp.expression_map(data3d, ROIs, markers)
    sp.expression_map(data3d, ROIs, markers, selected_types=['a', 'b', 'c'])
