import numpy as np
from spatialtis.abc import AnalysisBase
from spatialtis import Config

Config.verbose = False


def test_points_reader(data2d, data2d_keys, data3d, data3d_keys, data_wkt):
    ab = AnalysisBase(data=data2d, roi_key='ROI')
    assert np.array(ab.get_points()).shape[1] == 2

    ab = AnalysisBase(data=data2d_keys, roi_key='ROI', centroid_key=('x', 'y'))
    assert np.array(ab.get_points()).shape[1] == 2

    ab = AnalysisBase(data=data_wkt, roi_key='ROI')
    assert np.array(ab.get_points()).shape[1] == 2

    ab = AnalysisBase(data=data3d, roi_key='ROI')
    assert np.array(ab.get_points()).shape[1] == 3

    ab = AnalysisBase(data=data3d_keys, roi_key='ROI', centroid_key=('x', 'y', 'z'))
    assert np.array(ab.get_points()).shape[1] == 3


def test_exp_obs_naming(data2d):
    ab = AnalysisBase(data=data2d, roi_key="ROI")
    assert ab.exp_obs == ["ROI"]
    assert ab.roi_key == "ROI"

    ab = AnalysisBase(data=data2d, exp_obs=["ROI"])
    assert ab.exp_obs == ["ROI"]
    assert ab.roi_key == "ROI"
