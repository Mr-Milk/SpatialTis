import pytest

from spatialtis.preprocessing import read_ROIs

conditions = ["Patients", "Sample", "ROI"]


def test_read_rois(shared_datadir):
    data = read_ROIs((shared_datadir / 'Patients'), conditions, stacked=True)
    data.config_file((shared_datadir / 'metadata.csv'), channel_col="channels", marker_col="markers")
    data.to_anndata()


@pytest.mark.xfail
def test_concave(shared_datadir):
    data = read_ROIs((shared_datadir / 'Patients'), conditions, stacked=True)
    data.to_anndata(polygonize="concave", alpha=2.0)


def test_read_rois_mp(shared_datadir):
    data = read_ROIs((shared_datadir / 'Patients'), conditions, stacked=True)
    data.to_anndata(mp=True)
