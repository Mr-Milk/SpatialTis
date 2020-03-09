import pytest
from spatialtis.preprocessing import read_ROIs


DATA_ENTRY = "tests/data/Patients"
METADATA = "tests/data/metadata.csv"
conditions = ["Patients", "Sample", "ROI"]

data = read_ROIs(DATA_ENTRY, conditions, stacked=True)


def test_read_rois():
    data.config_file(METADATA, channel_col="channels", marker_col="markers")
    data.to_anndata()


@pytest.mark.xfail
def test_concave():
    data.to_anndata(polygonize="concave", alpha=2.0)


def test_read_rois_mp():
    data.to_anndata(mp=True)

