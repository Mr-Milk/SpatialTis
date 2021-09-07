"""
Setting Global config for whole processing level
"""
import warnings
from ast import literal_eval
from pathlib import Path
from typing import List, Optional, Sequence

import pandas as pd
from anndata import AnnData
from pyecharts.globals import WarningType, CurrentConfig, NotebookType
from rich.console import Console

WarningType.ShowWarning = False


class History(object):
    key = 'spatialtis_analysis_history'

    def __init__(self, data: AnnData):
        self.data = data
        if self.key in data.uns_keys():
            self.storage = data.uns[self.key]
        else:
            self.storage = {}
            data.uns[self.key] = self.storage

    def add(self, task_name: str, last_used_key: str):
        self.storage[task_name] = last_used_key
        self.data.uns[self.key] = self.storage


class _Config(object):
    def __init__(self):
        self._exp_obs: Optional[List[str]] = None
        self._cell_type_key: Optional[str] = None
        self._env: Optional[str] = None
        self._verbose: bool = True
        self.save_path: Optional[Path] = None
        self.progress_bar: bool = False

        self._roi_key: Optional[str] = None
        self.mp: bool = True

        # used key name to store info in anndata
        self.centroid_key: Optional[str] = None
        self.shape_key: Optional[str] = None
        self.marker_key: Optional[str] = None

    def __repr__(self):
        current_configs = [
            ['Multiprocessing', 'mp', str(self.mp)],
            ['Env', 'env', str(self.env)],
            ['Verbose', 'verbose', str(self.verbose)],
            ['Progress bar', 'progress_bar', str(self.progress_bar)],
            ["Auto save", "auto_save", str(self.save_path) if self.auto_save else str(self.auto_save)],
            ["Experiment observations", "exp_obs", str(self.exp_obs)],
            ["ROI key", "roi_key", str(self.roi_key)],
            ["Cell type key", "cell_type_key", str(self.cell_type_key)],
            ["Marker key", "marker_key", self.marker_key],
            ["Centroid key", "centroid_key", self.centroid_key],
            ["Shape key", "shape_key", self.shape_key],
        ]

        return pd.DataFrame(data=current_configs, columns=['Options', 'Attributes', 'Values'])

    def _to_dict(self):

        return dict(
            mp=self.mp,
            env=self.env,
            verbose=self.verbose,
            progress_bar=self.progress_bar,
            save_path=self.save_path,
            exp_obs=self.exp_obs,
            roi_key=self.roi_key,
            cell_type_key=self.cell_type_key,
            marker_key=self.marker_key,
            centroid_key=self.centroid_key,
            shape_key=self.shape_key,
        )

    @property
    def roi_key(self):
        return self._roi_key

    @roi_key.setter
    def roi_key(self, key):
        if self.exp_obs is None:
            self.exp_obs = [key]
        if key not in self.exp_obs:
            raise ValueError("The `roi_key` is not in your `exp_obs`")
        else:
            if self.exp_obs[-1] != key:
                exp_obs = self.exp_obs
                exp_obs.remove(key)
                exp_obs.append(key)
                self.exp_obs = exp_obs
            else:
                self._roi_key = key

    @property
    def exp_obs(self):
        return self._exp_obs

    @exp_obs.setter
    def exp_obs(self, obs):
        if isinstance(obs, (str, int, float)):
            self._exp_obs = [obs]
        elif isinstance(obs, Sequence):
            self._exp_obs = list(obs)
        elif obs is None:
            self._exp_obs = None
        else:
            raise TypeError(f"Couldn't set `exp_obs` with type {type(obs)}")
        self._roi_key = self._exp_obs[-1]

    @property
    def cell_type_key(self):
        return self._cell_type_key

    @cell_type_key.setter
    def cell_type_key(self, type_key):
        if isinstance(type_key, (str, int, float)):
            self._cell_type_key = type_key
        elif type_key is None:
            self._cell_type_key = None
        else:
            raise TypeError(f"Couldn't set CELL_TYPE_KEY with type {type(type_key)}")

    @property
    def env(self):
        return self._env

    @env.setter
    def env(self, env):
        self._env = env
        from bokeh.io import output_notebook

        output_notebook(hide_banner=True)

        if env is None:
            pass
        elif (env == "jupyter") | (env == "jupyter_lab"):
            CurrentConfig.NOTEBOOK_TYPE = NotebookType.JUPYTER_LAB
        elif env == "jupyter_notebook":
            CurrentConfig.NOTEBOOK_TYPE = NotebookType.JUPYTER_NOTEBOOK
        elif env == "nteract":
            CurrentConfig.NOTEBOOK_TYPE = NotebookType.NTERACT
        elif env == "zeppelin":
            CurrentConfig.NOTEBOOK_TYPE = NotebookType.ZEPPELIN
        elif env == "terminal":
            pass
        else:
            self._env = None
            warnings.warn("Unknown working environments, fallback to None", UserWarning)

    @property
    def verbose(self):
        return self._verbose

    @verbose.setter
    def verbose(self, v):
        if not isinstance(v, bool):
            raise ValueError("CONFIG.verbose only accept bool value")
        self._verbose = v
        self.progress_bar = v

    @property
    def auto_save(self):
        if self.save_path is not None:
            return True
        else:
            return False

    @auto_save.setter
    def auto_save(self, path):
        if isinstance(path, (str, Path)):
            path = Path(path)
            path.mkdir(exist_ok=True)
        elif isinstance(path, bool):
            if path:
                path = Path().cwd() / "spatialtis-result"
                path.mkdir(exist_ok=True)
            else:
                path = None
        else:
            raise TypeError("You can set a directory or set it True/False.")
        self.save_path = path

    def dumps(self, data: AnnData):
        data.uns['spatialtis_config'] = str(self._to_dict())

    def loads(self, data: AnnData):
        config = literal_eval(data.uns['spatialtis_config'])
        self.mp = config['mp']
        self.roi_key = config['roi_key']
        self.env = config['env']
        self.verbose = config['verbose']
        self.progress_bar = config['progress_bar']
        self.save_path = config['save_path']
        self.exp_obs = config['exp_obs']
        self.roi_key = config['roi_key']
        self.cell_type_key = config['cell_type_key']
        self.marker_key = config['marker_key']
        self.centroid_key = config['centroid_key']
        self.shape_key = config['shape_key']


# two state variable that don't need to reload
# need to init in the same place
Config = _Config()
console = Console()

if console.is_dumb_terminal:
    Config.env = None
elif console.is_jupyter:
    Config.env = "jupyter_lab"
elif console.is_terminal:
    Config.env = "terminal"
else:
    Config.env = None
    Config.progress_bar = False
