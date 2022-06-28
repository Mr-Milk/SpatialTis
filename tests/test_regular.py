import pytest

import spatialtis as st
import spatialtis.plotting as sp
from spatialtis import Config

Config.reset()

Config.verbose = False
Config.roi_key = "ROI"
Config.cell_type_key = "cell_type"
type_order = ["a", "b", "c"]


@pytest.mark.parametrize("method,kw", [
    ("id", dict(r=50)),
    ("morisita", dict(quad=(10, 10))),
    ("clark_evans", dict())
])
def test_cell_dispersion(data2d, method, kw):
    Config.loads(data2d)
    st.cell_dispersion(data2d, method=method, **kw)


def test_cell_dispersion_3d(data3d):
    Config.loads(data3d)
    st.cell_dispersion(data3d, r=50, method="id")


def test_cell_dispersion_plot(data2d):
    Config.loads(data2d)
    sp.cell_dispersion(data2d, type_order=type_order)
    sp.cell_dispersion(data2d, use="heatmap", type_order=type_order)


@pytest.mark.parametrize("method,kw", [
    ("leibovici", dict(d=10)),
    ("altieri", dict(cut=3)),
    ("shannon", dict())
])
def test_spatial_hetero(data2d, data3d, method, kw):
    Config.loads(data2d)
    st.spatial_heterogeneity(data2d, method=method, **kw)
    Config.loads(data3d)
    st.spatial_heterogeneity(data3d, method=method, **kw)


def test_spatial_hetero_plot(data2d):
    Config.loads(data2d)
    sp.spatial_heterogeneity(data2d)


def test_hotspot(data2d):
    Config.loads(data2d)
    st.hotspot(data2d, rect_side=(10, 10))

