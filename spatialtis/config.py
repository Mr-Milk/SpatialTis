"""
Setting Global config for whole processing level
"""
import platform
import sys
import warnings
from dataclasses import dataclass, field
from pathlib import Path
from typing import List, Optional, Sequence

from colorama import Fore
from pyecharts.globals import WarningType
from rich.table import Table

from spatialtis.console import console

WarningType.ShowWarning = False


@dataclass
class Analysis(object):
    task_name: str
    display_name: str
    export_key: str
    last_used_key: Optional[str] = field(default=None)


ANALYSIS = {
    "cell_components": Analysis(
        "cell_components", "Cell Components", "cell_components"
    ),
    "cell_density": Analysis("cell_density", "Cell Density", "cell_density"),
    "cell_co_occurrence": Analysis(
        "cell_co_occurrence", "Cell Co-occurrence", "co_occurrence"
    ),
    "cell_morphology": Analysis("cell_morphology", "Cell Morphology", "morphology"),
    "spatial_distribution": Analysis(
        "spatial_distribution", "Spatial Distribution", "spatial_distribution"
    ),
    "spatial_heterogeneity": Analysis(
        "spatial_heterogeneity", "Spatial Heterogeneity", "spatial_heterogeneity"
    ),
    "hotspot": Analysis("hotspot", "Hotspot Analysis", "hotspot"),
    "find_neighbors": Analysis(
        "find_neighbors", "Find Cell Neighbors", "cell_neighbors"
    ),
    "neighborhood_analysis": Analysis(
        "neighborhood_analysis", "Neighborhood Analysis", "neighborhood_analysis"
    ),
    "spatial_enrichment_analysis": Analysis(
        "spatial_enrichment_analysis",
        "Spatial Enrichment Analysis",
        "spatial_enrichment",
    ),
    "spatial_co_expression": Analysis(
        "spatial_co_expression", "Spatial Co-expression", "co_expression"
    ),
    "cell_community": Analysis("cell_community", "Cell Community", "community"),
    "NCDMarkers": Analysis(
        "NCDMarkers", "Neighbor cell dependent markers", "ncd_markers"
    ),
    "NMDMarkers": Analysis(
        "NMDMarkers", "Neighbor marker dependent markers", "nmd_markers"
    ),
    "prepare_svca": Analysis("prepare_scva", "Prepare files for running svca", "svca"),
}


class _CONFIG(object):
    def __init__(self):
        self._EXP_OBS: Optional[List[str]] = None
        self._CELL_TYPE_KEY: Optional[str] = None
        self._WORKING_ENV: Optional[str] = None
        self._VERBOSE: bool = True
        self.SAVE_PATH: Optional[Path] = None
        self.PBAR: bool = True

        self._ROI_KEY: Optional[str] = None
        self.OS: Optional[str] = None
        self.MP: bool = True

        # used key name to store info in anndata
        self.CENTROID_KEY: str = "centroid"
        self.AREA_KEY: str = "area"
        self.SHAPE_KEY: str = "cell_shape"
        self.ECCENTRICITY_KEY: str = "eccentricity"
        self.MARKER_KEY: str = "marker"
        self.NEIGHBORS_KEY: str = "cell_neighbors"
        self.pbar_format: str = f"{Fore.GREEN}{{desc}} {{bar}} {{percentage:3.0f}}% {{remaining}}|{{elapsed}}{Fore.RESET}"

    def __repr__(self):
        table = Table(title="Current configurations of SpatialTis")
        table.add_column("Option", style="cyan")
        table.add_column("Value", style="magenta")

        table.add_row("OS", self.OS)
        table.add_row("MULTI_PROCESSING", str(self.MP))
        table.add_row("WORKING_ENV", str(self.WORKING_ENV))
        table.add_row("VERBOSE", str(self.VERBOSE))
        table.add_row("PBAR", str(self.PBAR))
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
        from bokeh.io import output_notebook

        output_notebook(hide_banner=True)

        if env is None:
            pass
        elif env == "jupyter":
            pass
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
        self.PBAR = v

    @property
    def AUTO_SAVE(self):
        if self.SAVE_PATH is not None:
            return True
        else:
            return False

    @AUTO_SAVE.setter
    def AUTO_SAVE(self, path):
        if isinstance(path, (str, Path)):
            path = Path(path)
            path.mkdir(exist_ok=True)
        elif isinstance(path, bool):
            if path:
                path = Path().cwd() / "spatialti-result"
                path.mkdir(exist_ok=True)
            else:
                path = None
        else:
            raise TypeError("You can set a directory or set it True/False.")
        self.SAVE_PATH = path

    def pbar(self, **kwargs):
        pbar_config = dict(
            **kwargs,
            file=sys.stdout,
            disable=not self.PBAR,
            bar_format=self.pbar_format,
        )

        return pbar_config


CONFIG = _CONFIG()
# get system os
system_os = platform.system()
CONFIG.OS = system_os

if console.is_dumb_terminal:
    CONFIG.WORKING_ENV = None
    CONFIG.pbar_format = (
        f"{{desc}} {{bar}} {{percentage:3.0f}}% {{remaining}}|{{elapsed}}"
    )
elif console.is_jupyter:
    CONFIG.WORKING_ENV = "jupyter"
elif console.is_terminal:
    CONFIG.WORKING_ENV = "terminal"
else:
    CONFIG.WORKING_ENV = None
    CONFIG.PBAR = False
