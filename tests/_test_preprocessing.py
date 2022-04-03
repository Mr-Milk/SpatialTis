import sys

import pandas as pd
import pytest

from spatialtis import Config, read_ROIs

conditions = ["Patients", "Sample", "ROI"]
Config.env = None
Config.mp = False


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
    data.to_anndata(mp=False)
    data.to_anndata(polygonize="concave", alpha=1.0, mp=False)

    # when <py3.10, run this
    if sys.version_info[1] < 10:
        data.to_anndata(mp=True)
    else:
        return
