import pandas as pd
import pytest

from spatialtis import Config, read_ROIs

conditions = ["Patients", "Sample", "ROI"]
Config.env = None


def test_read_rois(tmpdir):
    var = pd.read_csv((tmpdir / "metadata.csv"))
    data = read_ROIs(
        (tmpdir / "Patients"),
        ["Patient", "Part", "ROI"],
        var,
        mask_pattern="mask",
        img_pattern="stacked",
    )
    pytest.data = data
    data.to_anndata()
    data.to_anndata(mp=True)


def test_concave(tmpdir):
    var = pd.read_csv((tmpdir / "metadata.csv"))
    data = read_ROIs(
        (tmpdir / "Patients" / "HPAP002" / "Body"),
        ["Part", "ROI"],
        var,
        mask_pattern="mask",
        img_pattern="stacked",
    )
    data.to_anndata(polygonize="concave", alpha=1.0)


@pytest.mark.xfail
def test_polygonzie_failed():
    data = pytest.data
    data.to_anndata(polygonize="con", alpha=1.0)
