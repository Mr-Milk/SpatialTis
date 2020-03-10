import pytest
from pathlib import Path

from spatialtis.preprocessing import read_ROIs

DATA_DIR = str(Path("tests/data/Patients").absolute())
METADATA_DIR = str(Path("tests/data/metadata.csv").absolute())

conditions = ["Patients", "Sample", "ROI"]

data = read_ROIs(DATA_DIR, conditions, stacked=True)


def test_read_rois():
    data.config_file(METADATA_DIR, channel_col="channels", marker_col="markers")
    data.to_anndata()


@pytest.mark.xfail
def test_concave():
    data.to_anndata(polygonize="concave", alpha=2.0)


def test_read_rois_mp():
    data.to_anndata(mp=True)
