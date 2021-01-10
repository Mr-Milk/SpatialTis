"""
Setting Global config for whole processing level
"""
import platform
import sys
import warnings
from typing import List, Optional, Sequence

from colorama import Fore
from pyecharts.globals import WarningType
from rich.table import Table

from spatialtis.console import console

WarningType.ShowWarning = False


class _CONFIG(object):
    def __init__(self):
        self._EXP_OBS: Optional[List[str]] = None
        self._CELL_TYPE_KEY: Optional[str] = None
        self._WORKING_ENV: Optional[str] = None
        self._VERBOSE: bool = True
        self.PBAR: bool = True

        self._ROI_KEY: Optional[str] = None
        self.OS: Optional[str] = None
        self.MULTI_PROCESSING: bool = False

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
        self.ncd_markers_key: str = "ncd_markers"
        self.nmd_markers_key: str = "nmd_markers"

    def __repr__(self):
        table = Table(title="Current configurations of SpatialTis")
        table.add_column("Option", style="cyan")
        table.add_column("Value", style="magenta")

        table.add_row("OS", self.OS)
        table.add_row("MULTI_PROCESSING", str(self.MULTI_PROCESSING))
        table.add_row("WORKING_ENV", str(self.WORKING_ENV))
        table.add_row("VERBOSE", str(self.VERBOSE))
        table.add_row("EXP_OBS", str(self.EXP_OBS))
        table.add_row("ROI_KEY", str(self.ROI_KEY))
        table.add_row("CELL_TYPE_KEY", str(self.CELL_TYPE_KEY))
        table.add_row("CENTROID_KEY", self.CENTROID_KEY)
        table.add_row("AREA_KEY", self.AREA_KEY)
        table.add_row("SHAPE_KEY", self.SHAPE_KEY)
        table.add_row("ECCENTRICITY_KEY", self.ECCENTRICITY_KEY)
        table.add_row("MARKER_KEY", self.MARKER_KEY)

        console.print(table)

        return ""

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
        self._WORKING_ENV = env
        from pyecharts.globals import CurrentConfig, NotebookType
        from bokeh.io import output_notebook

        output_notebook(hide_banner=True)

        if env is None:
            pass
        elif env == "jupyter_notebook":
            pass
        elif env == "jupyter_lab":
            CurrentConfig.NOTEBOOK_TYPE = NotebookType.JUPYTER_LAB
        elif env == "nteract":
            CurrentConfig.NOTEBOOK_TYPE = NotebookType.NTERACT
        elif env == "zeppelin":
            CurrentConfig.NOTEBOOK_TYPE = NotebookType.ZEPPELIN
        elif env == "terminal":
            pass
        else:
            self._WORKING_ENV = None
            warnings.warn("Unknown working environments, fallback to None", UserWarning)

    @property
    def VERBOSE(self):
        return self._VERBOSE

    @VERBOSE.setter
    def VERBOSE(self, v):
        if not isinstance(v, bool):
            raise ValueError("CONFIG.verbose only accept bool value")
        self._VERBOSE = v

    def pbar(self, **kwargs):
        pbar_config = dict(
            **kwargs,
            file=sys.stdout,
            disable=not self.PBAR,
            bar_format=f"{Fore.GREEN}{{desc}} {{bar}} {{percentage:3.0f}}% {{remaining}}|{{elapsed}}{Fore.RESET}",
        )

        return pbar_config


CONFIG = _CONFIG()
# get system os
system_os = platform.system()
CONFIG.OS = system_os

if console.is_dumb_terminal:
    CONFIG.WORKING_ENV = None
    CONFIG.PBAR = False
elif console.is_jupyter:
    CONFIG.WORKING_ENV = "jupyter_notebook"
elif console.is_terminal:
    CONFIG.WORKING_ENV = "terminal"
else:
    CONFIG.WORKING_ENV = None
    CONFIG.PBAR = False
