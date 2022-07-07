"""
Setting Global config for whole processing level
"""
import warnings
from anndata import AnnData
from ast import literal_eval
from pathlib import Path
from rich.console import Console
from rich.table import Table
from typing import List, Optional, Sequence

console = Console()


class History(object):
    key = "spatialtis_analysis_history"

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


class _Config:
    """Global configurations for spatialtis.

    Do not directly import this class, import the instance created for you.

    >>> from spatialtis import Config

    This allows you to set global configuration, so that you don't have to repeatly pass the same
    parameters in every function.

    To set a config, simple set the attribute.

    >>> Config.auto_save = True  # auto save to current directory
    >>> Config.auto_save = "my_spatialtis_result"  # save to custom directory

    To view your current configs.

    >>> Config.view()
           Current configurations of SpatialTis
    ┌─────────────────────────┬───────────────┬───────┐
    │ Options                 │ Attributes    │ Value │
    ├─────────────────────────┼───────────────┼───────┤
    │ Multiprocessing         │ mp            │ True  │
    │ Verbose                 │ verbose       │ True  │
    │ Progress bar            │ progress_bar  │ False │
    │ Auto save               │ auto_save     │ False │
    │ Experiment observations │ exp_obs       │ None  │
    │ ROI key                 │ roi_key       │       │
    │ Cell type key           │ cell_type_key │       │
    │ Marker key              │ marker_key    │       │
    │ Centroid key            │ centroid_key  │       │
    │ Shape key               │ shape_key     │       │
    └─────────────────────────┴───────────────┴───────┘


    Attributes
    ----------
    exp_obs : str or list of str, default: None
        **Required**. The columns in `.obs` that tells how your experiments organized,
        for example, you have different columns ['patients', 'sex', 'organ_part', 'roi_id'],
        the last one will be used as `roi_key`. You can override it by setting `roi_key`.
    centroid_key : str, default: None
        **Required**. The column in `.obs` or `.obsm` that store cell coordination,
        could be array-like or wkt format.
    roi_key : str, default: None
        Set the `roi_key`.
    shape_key : str, default: None
        The columns in `.obs` that store cell shape, must be in wkt format,
        use `spatialtis.wkt_shapes` to transform you data into wkt.
    cell_type_key : str, default: None
        The columns in `.obs` that store cell type name,
        some analyses require cell type name to proceed.
    marker_key : str, default: None
        The columns in `.var` that store protein/gene/transcript... name,
        if not specific, will use `.var.index`.
    mp : bool, default: True
        To turn on/off multiprocessing.
        From v0.5.0, this paramter has no effect.
    verbose : bool, default: True
        Control the printed message.
    progress_bar : bool, default True
        Control the progress bar.
    auto_save : bool or str or path, default: False
        Set to `True` will automatically save all images in
        `spatialtis_result` fold created at current directory,
        you can also pass a path to it.

    """

    def __init__(self):
        self._exp_obs: Optional[List[str]] = None
        self._cell_type_key: Optional[str] = None
        self._verbose: bool = True
        self.save_path: Optional[Path] = None
        self.progress_bar: bool = False

        self._roi_key: Optional[str] = None
        self.mp: bool = True

        # used key name to store info in anndata
        self._centroid_key: Optional[str] = None
        self.shape_key: Optional[str] = None
        self.marker_key: Optional[str] = None

    def view(self):
        """Print a table with current configurations"""
        table = Table(title="Current configurations of SpatialTis")
        table.add_column("Options", style="bright_black")
        table.add_column("Attributes", style="cyan")
        table.add_column("Value", style="magenta")
        # add content
        table.add_row("Multiprocessing", "mp", str(self.mp))
        table.add_row("Verbose", "verbose", str(self.verbose))
        table.add_row("Progress bar", "progress_bar", str(self.progress_bar))
        table.add_row(
            "Auto save",
            "auto_save",
            str(self.save_path) if self.auto_save else str(self.auto_save),
        )
        render_exp_obs = "" if self.exp_obs is None else str(self.exp_obs)
        table.add_row("Experiment observations", "exp_obs", render_exp_obs)
        table.add_row("ROI key", "roi_key", self.roi_key)
        table.add_row("Cell type key", "cell_type_key", self.cell_type_key)
        table.add_row("Marker key", "marker_key", self.marker_key)
        table.add_row("Centroid key", "centroid_key", self.centroid_key)
        table.add_row("Shape key", "shape_key", self.shape_key)
        console.print(table)

    def __repr__(self):
        self.view()
        return "SpatialTis Global Configurations"

    def _to_dict(self):

        return dict(
            mp=self.mp,
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
        self._roi_key = key
        if self.exp_obs is None:
            if key is not None:
                self.exp_obs = [key]
        if key not in self.exp_obs:
            self.exp_obs = [*self._exp_obs, key]
        else:
            if self.exp_obs[-1] != key:
                exp_obs = self.exp_obs
                exp_obs.remove(key)
                exp_obs.append(key)
                self.exp_obs = exp_obs

    @property
    def exp_obs(self):
        return self._exp_obs

    @exp_obs.setter
    def exp_obs(self, obs):
        if isinstance(obs, (str, int, float)):
            self._exp_obs = [obs]
            self._roi_key = self._exp_obs[-1]
        elif isinstance(obs, Sequence):
            if len(obs) == 0:
                self._exp_obs = None
            else:
                self._exp_obs = list(obs)
                self._roi_key = self._exp_obs[-1]
        elif obs is None:
            self._exp_obs = None
        else:
            raise TypeError(f"Couldn't set `exp_obs` with type {type(obs)}")

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
            raise TypeError(f"Couldn't set cell_type_key with type {type(type_key)}")

    @property
    def centroid_key(self):
        return self._centroid_key

    @centroid_key.setter
    def centroid_key(self, key):
        if isinstance(key, (str, int, float)):
            self._centroid_key = key
        elif isinstance(key, Sequence):
            if len(key) == 0:
                self._centroid_key = None
            else:
                self._centroid_key = key
        elif key is None:
            self._centroid_key = None
        else:
            raise TypeError(f"Couldn't set centroid_key with type {type(key)}")

    @property
    def verbose(self):
        return self._verbose

    @verbose.setter
    def verbose(self, v):
        if not isinstance(v, bool):
            raise ValueError("Config.verbose only accept bool value")
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
        """Save configurations to anndata

        Parameters
        ----------
        data : AnnData
            The `AnnData` object to save the Config

        """
        data.uns["spatialtis_config"] = str(self._to_dict())

    def loads(self, data: AnnData):
        """Load configurations from anndata

        Parameters
        ----------
        data : AnnData
            The `AnnData` object to load the Config

        """
        store_key = "spatialtis_config"
        if store_key in data.uns_keys():
            config = literal_eval(data.uns[store_key])
            self.reset()
            self.mp = config["mp"]
            self.exp_obs = config["exp_obs"]
            self.roi_key = config["roi_key"]
            self.verbose = config["verbose"]
            self.progress_bar = config["progress_bar"]
            self.save_path = config["save_path"]
            self.roi_key = config["roi_key"]
            self.cell_type_key = config["cell_type_key"]
            self.marker_key = config["marker_key"]
            self.centroid_key = config["centroid_key"]
            self.shape_key = config["shape_key"]
        else:
            warnings.warn("No config is found")

    def reset(self):
        """Reset to default"""

        self._verbose = True
        self._exp_obs = None
        self._roi_key = None
        self._cell_type_key = None

        self.mp = True
        self.progress_bar = True
        self.auto_save = False
        self.marker_key = None
        self.centroid_key = None
        self.shape_key = None


# two state variable that don't need to reload
# need to init in the same place
Config = _Config()
