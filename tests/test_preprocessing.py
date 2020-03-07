import pytest
from spatialtis.preprocessing import read_ROIs

DATA_ENTRY = "../data/shrink_data/Patients"
METADATA = "../data/shrink_data/metadata.csv"

conditions = ["Patients", "Sample", "ROI"]

data = read_ROIs(DATA_ENTRY, conditions, stacked=True)

data.config_file(METADATA, channel_col="channels", marker_col="markers")

data.to_anndata()


@pytest.mark.xfail
def test_concave():
    data.to_anndata(polygonize="concave", alpha=2.0)


data.to_anndata(mp=True)
