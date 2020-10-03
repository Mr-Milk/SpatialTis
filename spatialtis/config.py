"""
Setting Global config for whole processing level
"""
import platform
import sys
import warnings
from typing import Any, Mapping, Optional, Sequence

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
        self._EXP_OBS: Optional[Sequence[str]] = None
        self._CELL_TYPE_KEY: Optional[str] = None
        self._WORKING_ENV: Optional[str] = None
        self._CPU_ALLOC: Optional[int] = None
        self._CONFIG_VERBOSE: Any = _VERBOSE()

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
            f"CPU_ALLOC: {self.CPU_ALLOC}",
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
    def CPU_ALLOC(self):
        return self._CPU_ALLOC

    @CPU_ALLOC.setter
    def CPU_ALLOC(self, num):
        # auto run ray.init when running on Linux and MacOS
        if not isinstance(num, int):
            raise ValueError("The number of CPU must be integer.")

        self._CPU_ALLOC = num

        import logging
        import ray

        ray.init(
            logging_level=logging.FATAL,
            ignore_reinit_error=True,
            num_cpus=self._CPU_ALLOC,
        )

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

try:
    import logging
    import ray

    ray.init(
        logging_level=logging.FATAL, ignore_reinit_error=True,
    )
except Exception:
    warnings.warn("Initiate ray failed, parallel processing not available.")

# can be override in plotting function
# ['jupyter', 'zepplin', None]

# ==============SOME INFO==================
ISOTOPES_NAME: Sequence = [
    "H",
    "He",
    "Li",
    "Be",
    "B",
    "C",
    "N",
    "O",
    "F",
    "Ne",
    "Na",
    "Mg",
    "Al",
    "Si",
    "P",
    "S",
    "Cl",
    "Ar",
    "K",
    "Ca",
    "Sc",
    "Ti",
    "V",
    "Cr",
    "Mn",
    "Fe",
    "Co",
    "Ni",
    "Cu",
    "Zn",
    "Ga",
    "Ge",
    "As",
    "Se",
    "Br",
    "Kr",
    "Rb",
    "Sr",
    "Y",
    "Zr",
    "Nb",
    "Mo",
    "Tc",
    "Ru",
    "Rh",
    "Pd",
    "Ag",
    "Cd",
    "In",
    "Sn",
    "Sb",
    "Te",
    "I",
    "Xe",
    "Cs",
    "Ba",
    "La",
    "Ce",
    "Pr",
    "Nd",
    "Pm",
    "Sm",
    "Eu",
    "Gd",
    "Tb",
    "Dy",
    "Ho",
    "Er",
    "Tm",
    "Yb",
    "Lu",
    "Hf",
    "Ta",
    "W",
    "Re",
    "Os",
    "Ir",
    "Pt",
    "Au",
    "Hg",
    "Tl",
    "Pb",
    "Bi",
    "Po",
    "At",
    "Rn",
    "Fr",
    "Ra",
    "Ac",
    "Th",
    "Pa",
    "U",
    "Np",
    "Pu",
    "Am",
    "Cm",
    "Bk",
    "Cf",
    "Es",
    "Fm",
    "Md",
    "No",
    "Lr",
    "Rf",
    "Db",
    "Sg",
    "Bh",
    "Hs",
    "Mt",
    "Ds",
    "Rg",
    "Cn",
    "Nh",
    "Fl",
    "Mc",
    "Lv",
    "Ts",
    "Og",
]

ISOTOPES_MASS_NUMBER_MAP: Mapping = {
    "H": [1, 2, 3],
    "He": [3, 4],
    "Li": [6, 7],
    "Be": [9],
    "B": [10, 11],
    "C": [12, 13, 14],
    "N": [14, 15],
    "O": [16, 17, 18],
    "F": [19],
    "Ne": [20, 21, 22],
    "Na": [23, 22],
    "Mg": [24, 25, 26],
    "Al": [27],
    "Si": [28, 29, 30],
    "P": [31],
    "S": [32, 33, 34, 36],
    "Cl": [35, 37, 36],
    "Ar": [36, 38, 40, 39],
    "K": [39, 40, 41],
    "Ca": [40, 42, 43, 44, 46, 48, 41],
    "Sc": [45],
    "Ti": [46, 47, 48, 49, 50],
    "V": [50, 51],
    "Cr": [50, 52, 53, 54],
    "Mn": [55, 53],
    "Fe": [54, 56, 57, 58],
    "Co": [59, 60],
    "Ni": [58, 60, 61, 62, 64],
    "Cu": [63, 65],
    "Zn": [64, 66, 67, 68, 70],
    "Ga": [69, 71],
    "Ge": [70, 72, 73, 74, 76],
    "As": [75],
    "Se": [74, 76, 77, 78, 80, 82, 79],
    "Br": [79, 81],
    "Kr": [78, 80, 82, 83, 84, 86, 85],
    "Rb": [85, 87],
    "Sr": [84, 86, 87, 88],
    "Y": [89],
    "Zr": [90, 91, 92, 94, 96],
    "Nb": [93],
    "Mo": [92, 94, 95, 96, 97, 98, 100],
    "Tc": [98, 97, 99],
    "Ru": [96, 98, 99, 100, 101, 102, 104],
    "Rh": [103],
    "Pd": [102, 104, 105, 106, 108, 110],
    "Ag": [107, 109],
    "Cd": [106, 108, 110, 111, 112, 113, 114, 116],
    "In": [113, 115],
    "Sn": [112, 114, 115, 116, 117, 118, 119, 120, 122, 124],
    "Sb": [121, 123, 125],
    "Te": [120, 122, 123, 124, 125, 126, 128, 130],
    "I": [127, 129],
    "Xe": [124, 126, 128, 129, 130, 131, 132, 134, 136],
    "Cs": [133, 134, 135, 137],
    "Ba": [130, 132, 134, 135, 136, 137, 138, 133],
    "La": [138, 139, 137],
    "Ce": [136, 138, 140, 142],
    "Pr": [141],
    "Nd": [142, 143, 144, 145, 146, 148, 150],
    "Pm": [145, 146, 147],
    "Sm": [144, 147, 148, 149, 150, 152, 154, 151],
    "Eu": [151, 153, 152, 154, 155],
    "Gd": [152, 154, 155, 156, 157, 158, 160],
    "Tb": [159, 157, 160],
    "Dy": [156, 158, 160, 161, 162, 163, 164],
    "Ho": [165],
    "Er": [162, 164, 166, 167, 168, 170],
    "Tm": [169, 171],
    "Yb": [168, 170, 171, 172, 173, 174, 176],
    "Lu": [175, 176, 173, 174],
    "Hf": [174, 176, 177, 178, 179, 180],
    "Ta": [180, 181],
    "W": [180, 182, 183, 184, 186],
    "Re": [185, 187],
    "Os": [184, 186, 187, 188, 189, 190, 192],
    "Ir": [191, 193],
    "Pt": [190, 192, 194, 195, 196, 198],
    "Au": [197],
    "Hg": [196, 198, 199, 200, 201, 202, 204],
    "Tl": [203, 205, 204],
    "Pb": [204, 206, 207, 208],
    "Bi": [209, 207],
    "Po": [208, 209, 210],
    "At": [210, 211],
    "Rn": [210, 211, 222],
    "Fr": [212, 222, 223],
    "Ra": [226, 228],
    "Ac": [225, 227],
    "Th": [230, 232, 229],
    "Pa": [231, 233],
    "U": [233, 234, 235, 238, 236],
    "Np": [236, 237],
    "Pu": [238, 239, 240, 241, 242, 244],
    "Am": [241, 243],
    "Cm": [243, 244, 245, 246, 247, 248],
    "Bk": [247, 249],
    "Cf": [249, 250, 251, 252],
    "Es": [252, 254],
    "Fm": [253, 257],
    "Md": [258, 260],
    "No": [255, 259],
    "Lr": [261, 262],
    "Rf": [265, 267],
    "Db": [268, 270],
    "Sg": [269, 271],
    "Bh": [270, 274],
    "Hs": [269, 270],
    "Mt": [276, 278],
    "Ds": [280, 281],
    "Rg": [281, 282],
    "Cn": [283, 285],
    "Nh": [285, 286],
    "Fl": [287, 288, 289],
    "Mc": [288, 289, 290],
    "Lv": [291, 292, 293],
    "Ts": [293, 294],
    "Og": [294],
}
