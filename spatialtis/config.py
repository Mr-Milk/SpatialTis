"""
Setting Global config for whole processing level
"""
import platform
import sys
import warnings
from typing import List, Optional, Sequence

from colorama import Fore
from pyecharts.globals import WarningType

WarningType.ShowWarning = False


class _VERBOSE(object):
    def __init__(self):
        self.PBAR = True
        self.ANNDATA = True
        self.INFO = True

    def __repr__(self):
        info = [f"PBAR: {self.PBAR}", f"ANNDATA: {self.ANNDATA}", f"INFO: {self.INFO}"]
        return "\n".join(info)


class _CONFIG(object):
    def __init__(self):
        self._EXP_OBS: Optional[List[str]] = None
        self._CELL_TYPE_KEY: Optional[str] = None
        self._WORKING_ENV: Optional[str] = None
        self._CPU_ALLOC: Optional[int] = None
        self._CONFIG_VERBOSE: _VERBOSE = _VERBOSE()

        self._ROI_KEY: Optional[str] = None
        self.OS: Optional[str] = None
        self.MULTI_PROCESSING: bool = False

        # set tqdm bar foramt
        self._INBUILT_PBAR_FORMAT = "%s{l_bar}%s{bar}%s{r_bar}%s" % (
            Fore.GREEN,
            Fore.CYAN,
            Fore.GREEN,
            Fore.RESET,
        )
        self.PBAR_FORMAT = self._INBUILT_PBAR_FORMAT

        # used key name to store info in anndata
        self.CENTROID_KEY: str = "centroid"
        self.AREA_KEY: str = "area"
        self.SHAPE_KEY: str = "cell_shape"
        self.ECCENTRICITY_KEY: str = "eccentricity"
        self.MARKER_KEY: str = "marker"

        # export key, the key name used to store results, private to user
        # statistic part
        self.cell_components_key: str = "cell_components"
        self.cell_co_occurrence_key: str = "cell_co_occurrence"
        self.cell_density_key: str = "cell_density"
        self.cell_morphology_key: str = "cell_morphology"

        # plotting part
        self.spatial_distribution_key: str = "spatial_distribution"
        self.spatial_heterogeneity_key: str = "spatial_heterogeneity"
        self.hotspot_key: str = "hotspot"
        self.community_key: str = "communities"
        self.neighborhood_analysis_key: str = "neighborhood_analysis"
        self.spatial_enrichment_analysis_key: str = "spatial_enrichment_analysis"
        self.spatial_enrichment_analysis_layers_key: str = "markers_sign"
        self.neighbors_key: str = "cell_neighbors"
        self.neighbors_count_key: str = "neighbors_count"
        self.exp_neighcell_key: str = "exp_neighcell"
        self.exp_neighexp_key: str = "exp_neighexp"

    def __repr__(self):
        head = "Current configurations of SpatialTis:"
        info = [
            "\n======CONFIGURATIONS======",
            f"EXP_OBS: {self.EXP_OBS}",
            f"OS: {self.OS}",
            f"MULTI_PROCESSING: {self.MULTI_PROCESSING}",
            f"WORKING_ENV: {self.WORKING_ENV}",
            f"VERBOSE.PBAR: {self.VERBOSE.PBAR}",
            f"VERBOSE.ANNDATA: {self.VERBOSE.ANNDATA}",
            f"VERBOSE.INFO: {self.VERBOSE.INFO}",
            "\n======KEYS======",
            f"CELL_TYPE_KEY: {self.CELL_TYPE_KEY}",
            f"CENTROID_KEY: {self.CENTROID_KEY}",
            f"AREA_KEY: {self.AREA_KEY}",
            f"SHAPE_KEY: {self.SHAPE_KEY}",
            f"ECCENTRICITY_KEY: {self.ECCENTRICITY_KEY}",
            f"MARKER_KEY: {self.MARKER_KEY}",
        ]
        info = "\n".join(info)
        return f"{head}\n{info}"

    def tqdm(self, **kwargs):
        all_kwargs = dict(
            unit="ROI",
            bar_format=self.PBAR_FORMAT,
            disable=(not self.VERBOSE.PBAR),
            file=sys.stdout,
        )
        for k, v in kwargs.items():
            all_kwargs[k] = v
        return all_kwargs

    @property
    def ROI_KEY(self):
        return self._ROI_KEY

    @ROI_KEY.setter
    def ROI_KEY(self, key):
        if self.EXP_OBS is None:
            raise ValueError("Please set the EXP_OBS first")
        if key not in self.EXP_OBS:
            raise ValueError("The ROI_KEY is not in your EXP_OBS")
        else:
            if self.EXP_OBS[-1] != key:
                exp_obs = self.EXP_OBS
                exp_obs.remove(key)
                exp_obs.append(key)
                self.EXP_OBS = exp_obs
            else:
                self._ROI_KEY = key

    @property
    def EXP_OBS(self):
        return self._EXP_OBS

    @EXP_OBS.setter
    def EXP_OBS(self, obs):
        if isinstance(obs, (str, int, float)):
            self._EXP_OBS = [obs]
            self._ROI_KEY = self._EXP_OBS[-1]
        elif isinstance(obs, Sequence):
            self._EXP_OBS = list(obs)
            self._ROI_KEY = self._EXP_OBS[-1]
        elif obs is None:
            self._EXP_OBS = None
        else:
            raise TypeError(f"Couldn't set EXP_OBS with type {type(obs)}")

    @property
    def CELL_TYPE_KEY(self):
        return self._CELL_TYPE_KEY

    @CELL_TYPE_KEY.setter
    def CELL_TYPE_KEY(self, type_key):
        if isinstance(type_key, (str, int, float)):
            self._CELL_TYPE_KEY = type_key
        elif type_key is None:
            self._CELL_TYPE_KEY = None
        else:
            raise TypeError(f"Couldn't set CELL_TYPE_KEY with type {type(type_key)}")

    @property
    def WORKING_ENV(self):
        return self._WORKING_ENV

    @WORKING_ENV.setter
    def WORKING_ENV(self, env):
        if env not in [
            "jupyter_notebook",
            "jupyter_lab",
            "nteract",
            "zeppelin",
            "interactive",
            None,
        ]:
            warnings.warn("Unknown working environments", UserWarning)
        self._WORKING_ENV = env
        if env is None:
            self.VERBOSE.PBAR = False
        else:
            self.VERBOSE.PBAR = True
            self.PBAR_FORMAT = self._INBUILT_PBAR_FORMAT

            # TODO: Currently we have to turn of tqdm when run in notebook env in Windows,
            if (CONFIG.OS == "Windows") & (env != "interactive"):
                CONFIG.VERBOSE.PBAR = False

            from pyecharts.globals import CurrentConfig, NotebookType
            from bokeh.io import output_notebook

            output_notebook(hide_banner=True)

            if env == "jupyter_notebook":
                pass
            elif env == "jupyter_lab":
                CurrentConfig.NOTEBOOK_TYPE = NotebookType.JUPYTER_LAB
            elif env == "nteract":
                CurrentConfig.NOTEBOOK_TYPE = NotebookType.NTERACT
            elif env == "zeppelin":
                CurrentConfig.NOTEBOOK_TYPE = NotebookType.ZEPPELIN

    @property
    def VERBOSE(self):
        return self._CONFIG_VERBOSE

    @VERBOSE.setter
    def VERBOSE(self, v):
        if not isinstance(v, bool):
            raise ValueError("CONFIG.verbose only accept bool value")

        if v:
            self._CONFIG_VERBOSE.PBAR = True
            self._CONFIG_VERBOSE.ANNDATA = True
            self._CONFIG_VERBOSE.INFO = True
        else:
            self._CONFIG_VERBOSE.PBAR = False
            self._CONFIG_VERBOSE.ANNDATA = False
            self._CONFIG_VERBOSE.INFO = False


CONFIG = _CONFIG()
# get system os
system_os = platform.system()
CONFIG.OS = system_os
# I have no idea how to fix tqdm multi lines on jupyter environment for windows
if CONFIG.OS == "windows":
    CONFIG.VERBOSE.PBAR = False
