import numpy as np
from spatialtis.abc import AnalysisBase
from spatialtis import Config


def test_points_reader(data2d, data2d_keys, data3d, data3d_keys, data_wkt):
    Config.loads(data2d)
    ab = AnalysisBase(data=data2d)
    assert np.array(ab.get_centroids()).shape[1] == 2

    Config.loads(data2d_keys)
    ab = AnalysisBase(data=data2d_keys)
    assert np.array(ab.get_centroids()).shape[1] == 2

    Config.loads(data_wkt)
    ab = AnalysisBase(data=data_wkt)
    assert np.array(ab.get_centroids()).shape[1] == 2

    Config.loads(data3d)
    ab = AnalysisBase(data=data3d)
    assert np.array(ab.get_centroids()).shape[1] == 3

    Config.loads(data3d_keys)
    ab = AnalysisBase(data=data3d_keys)
    assert np.array(ab.get_centroids()).shape[1] == 3


def test_exp_obs_naming(data2d):
    ab = AnalysisBase(data=data2d, roi_key="ROI")
    assert ab.exp_obs == ["ROI"]
    assert ab.roi_key == "ROI"

    ab = AnalysisBase(data=data2d, exp_obs=["ROI"])
    assert ab.exp_obs == ["ROI"]
    assert ab.roi_key == "ROI"
