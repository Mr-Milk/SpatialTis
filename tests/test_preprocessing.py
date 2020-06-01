import pandas as pd
import pytest

from spatialtis import read_ROIs

conditions = ["Patients", "Sample", "ROI"]


def test_read_rois(shared_datadir):
    var = pd.read_csv((shared_datadir / 'metadata.csv'))
    data = read_ROIs((shared_datadir / 'Patients'), ['Patient', 'Part', 'ROI'], var,
                     mask_pattern="mask", img_pattern="stacked")
    pytest.data = data
    data.to_anndata()


def test_read_rois_mp():
    data = pytest.data
    data.to_anndata(mp=True)


def test_concave():
    data = pytest.data
    data.to_anndata(polygonize="concave", alpha=1.0, mp=True)
