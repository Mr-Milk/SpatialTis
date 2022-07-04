import warnings
from anndata import AnnData

from spatialtis import Config
from spatialtis.abc import AnalysisBase
from spatialtis.utils import default_args


def list_roi(data: AnnData, **kwargs):
    ab = AnalysisBase(data, verbose=False, **kwargs)
    return [roi_name for [roi_name] in ab.iter_roi(disable_pbar=True)]


def make_roi_unique(data: AnnData, write_config=True, key=None, **kwargs):
    ab = AnalysisBase(data, verbose=False, **kwargs)
    if ab.is_rois_name_unique(warn=False):
        msg = "Your ROI names are already unique."
        warnings.warn(msg)
    roi_names = []
    for roi_name, ix in ab.iter_roi(fields=["index"], disable_pbar=True):
        use_name = "_".join(roi_name)
        roi_names += [use_name for _ in range(len(ix))]
    key = default_args(key, f"{ab.roi_key}_unique")
    Config.exp_obs = [*Config.exp_obs[:-1], key]
    data.obs[key] = roi_names
