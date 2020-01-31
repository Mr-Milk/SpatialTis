import re
import numpy as np
import pandas as pd
import anndata as ad
from skimage.io import imread
from pathlib import Path
from skimage.external import tifffile
from collections import OrderedDict

from typing import Iterable, List, Union

from .geom import polygonize_cells

ISOTOPES_NAME = ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar',
                 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br',
                 'Kr', 'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te',
                 'I', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm',
                 'Yb', 'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn',
                 'Fr', 'Ra', 'Ac', 'Th', 'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr',
                 'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt', 'Ds', 'Rg', 'Cn', 'Nh', 'Fl', 'Mc', 'Lv', 'Ts', 'Og']
ISOTOPES_MASS_NUMBER_MAP = {'H': [1, 2, 3],
                            'He': [3, 4],
                            'Li': [6, 7],
                            'Be': [9],
                            'B': [10, 11],
                            'C': [12, 13, 14],
                            'N': [14, 15],
                            'O': [16, 17, 18],
                            'F': [19],
                            'Ne': [20, 21, 22],
                            'Na': [23, 22],
                            'Mg': [24, 25, 26],
                            'Al': [27],
                            'Si': [28, 29, 30],
                            'P': [31],
                            'S': [32, 33, 34, 36],
                            'Cl': [35, 37, 36],
                            'Ar': [36, 38, 40, 39],
                            'K': [39, 40, 41],
                            'Ca': [40, 42, 43, 44, 46, 48, 41],
                            'Sc': [45],
                            'Ti': [46, 47, 48, 49, 50],
                            'V': [50, 51],
                            'Cr': [50, 52, 53, 54],
                            'Mn': [55, 53],
                            'Fe': [54, 56, 57, 58],
                            'Co': [59, 60],
                            'Ni': [58, 60, 61, 62, 64],
                            'Cu': [63, 65],
                            'Zn': [64, 66, 67, 68, 70],
                            'Ga': [69, 71],
                            'Ge': [70, 72, 73, 74, 76],
                            'As': [75],
                            'Se': [74, 76, 77, 78, 80, 82, 79],
                            'Br': [79, 81],
                            'Kr': [78, 80, 82, 83, 84, 86, 85],
                            'Rb': [85, 87],
                            'Sr': [84, 86, 87, 88],
                            'Y': [89],
                            'Zr': [90, 91, 92, 94, 96],
                            'Nb': [93],
                            'Mo': [92, 94, 95, 96, 97, 98, 100],
                            'Tc': [98, 97, 99],
                            'Ru': [96, 98, 99, 100, 101, 102, 104],
                            'Rh': [103],
                            'Pd': [102, 104, 105, 106, 108, 110],
                            'Ag': [107, 109],
                            'Cd': [106, 108, 110, 111, 112, 113, 114, 116],
                            'In': [113, 115],
                            'Sn': [112, 114, 115, 116, 117, 118, 119, 120, 122, 124],
                            'Sb': [121, 123, 125],
                            'Te': [120, 122, 123, 124, 125, 126, 128, 130],
                            'I': [127, 129],
                            'Xe': [124, 126, 128, 129, 130, 131, 132, 134, 136],
                            'Cs': [133, 134, 135, 137],
                            'Ba': [130, 132, 134, 135, 136, 137, 138, 133],
                            'La': [138, 139, 137],
                            'Ce': [136, 138, 140, 142],
                            'Pr': [141],
                            'Nd': [142, 143, 144, 145, 146, 148, 150],
                            'Pm': [145, 146, 147],
                            'Sm': [144, 147, 148, 149, 150, 152, 154, 151],
                            'Eu': [151, 153, 152, 154, 155],
                            'Gd': [152, 154, 155, 156, 157, 158, 160],
                            'Tb': [159, 157, 160],
                            'Dy': [156, 158, 160, 161, 162, 163, 164],
                            'Ho': [165],
                            'Er': [162, 164, 166, 167, 168, 170],
                            'Tm': [169, 171],
                            'Yb': [168, 170, 171, 172, 173, 174, 176],
                            'Lu': [175, 176, 173, 174],
                            'Hf': [174, 176, 177, 178, 179, 180],
                            'Ta': [180, 181],
                            'W': [180, 182, 183, 184, 186],
                            'Re': [185, 187],
                            'Os': [184, 186, 187, 188, 189, 190, 192],
                            'Ir': [191, 193],
                            'Pt': [190, 192, 194, 195, 196, 198],
                            'Au': [197],
                            'Hg': [196, 198, 199, 200, 201, 202, 204],
                            'Tl': [203, 205, 204],
                            'Pb': [204, 206, 207, 208],
                            'Bi': [209, 207],
                            'Po': [208, 209, 210],
                            'At': [210, 211],
                            'Rn': [210, 211, 222],
                            'Fr': [212, 222, 223],
                            'Ra': [226, 228],
                            'Ac': [225, 227],
                            'Th': [230, 232, 229],
                            'Pa': [231, 233],
                            'U': [233, 234, 235, 238, 236],
                            'Np': [236, 237],
                            'Pu': [238, 239, 240, 241, 242, 244],
                            'Am': [241, 243],
                            'Cm': [243, 244, 245, 246, 247, 248],
                            'Bk': [247, 249],
                            'Cf': [249, 250, 251, 252],
                            'Es': [252, 254],
                            'Fm': [253, 257],
                            'Md': [258, 260],
                            'No': [255, 259],
                            'Lr': [261, 262],
                            'Rf': [265, 267],
                            'Db': [268, 270],
                            'Sg': [269, 271],
                            'Bh': [270, 274],
                            'Hs': [269, 270],
                            'Mt': [276, 278],
                            'Ds': [280, 281],
                            'Rg': [281, 282],
                            'Cn': [283, 285],
                            'Nh': [285, 286],
                            'Fl': [287, 288, 289],
                            'Mc': [288, 289, 290],
                            'Lv': [291, 292, 293],
                            'Ts': [293, 294],
                            'Og': [294]}


def mask2cells(
        mask_img: Union[Path, str],
        ignore_bg: bool = True
) -> Iterable[List]:
    """
    Parameters
    mask_img: Path/str, the path to mask img
    ignore_bg: bool, ignore background in mask
    """
    # read mask image
    mask = imread(mask_img)
    # number of cells
    counts = np.unique(mask)
    # create list for each cell
    cells = [[] for i in range(0, len(counts))]
    # find the exact points that belong to each cell
    it = np.nditer(mask, flags=['multi_index'])
    while not it.finished:
        cells[int(it[0])].append(it.multi_index)
        it.iternext()
    # usually 0 for background
    if ignore_bg:
        return cells[1:]
    else:
        return cells


def get_cell_exp_single_channel(
        channel: Iterable,
        cells: Iterable,
        method: str = 'mean'
) -> list:
    """
    Parameters
    channel: list or ndarray, matrix info of the channel
    cells: list or ndarray, each element contains points for each cell
    method: str, ('mean' / 'median' / 'sum' / ...)  any numpy method to compute expression level of single cell
    """
    cells_density = list()
    for cell in cells:
        density = channel[[i[0] for i in cell], [i[1] for i in cell]]
        exec(f'cells_density.append(np.{method}(density))')
    return cells_density


def get_cell_exp_stack(
        stack: Iterable,
        cells: Iterable,
        method: str = 'mean'
) -> np.ndarray:
    """
    Parameters
    channel: list or ndarray, matrix info of the channel
    cells: list or ndarray, each element contains points for each cell
    method: str, ('mean' / 'median' / 'sum' / ...)  any numpy method to compute expression level of single cell
    """
    cells_density = list()
    for cell in cells:
        cell_pixels = stack[:, [i[0] for i in cell], [i[1] for i in cell]]
        exec(f'cells_density.append([np.{method}(density) for density in cell_pixels])')
    # each secondary array as a channel
    return np.array(cells_density)


class read_all_ROIs:

    def __init__(self, entry, conditions, mask_pattern="*mask*"):
        """
        Organize your folder as your experiment conditions.
        We will exhaust the directory.
        Length of condition_name must be the same as the depth of your folders.
    
        Example:

        condition_name = ['Patient', 'Sample', 'ROI']

        Entry
        |------Patient1
            |------Sample1
                |------ROI1
                |------ROI2
            |------Sample2
                |------ROI1
                |------ROI2
            |------Sample3
                |------ROI1
                |------ROI2
        |------Patient2
            |------Sample1
            |------Sample2
            |------Sample3
        |------Patient3
            |------Sample1
            |------Sample2
            |------Sample3
        
        """

        self._mask_pattern = mask_pattern
        self._channels = list()
        self._markers = dict()
        self._channels_files = dict()

        self._tree = []
        self._obs_name = conditions
        self._depth = len(conditions)
        self._exhaust_dir(entry)
        self._obs = [p.parts[-self._depth:] for p in self._tree]
        self._var = None

    # walk through the directory, until there is no directory
    def _exhaust_dir(self, path):
        d = [f for f in Path(path).iterdir() if f.is_dir()]
        for f in d:
            self._tree.append(f)
            if f.parent in self._tree:
                self._tree.remove(f.parent)
            self._exhaust_dir(f)

    @property
    def tree(self):
        return self._obs

    @property
    def vars(self):
        return self._var

    def config_file(self, metadata, channel_col=None, marker_col=None, sep=','):
        meta = pd.read_csv(metadata, sep=sep)
        channels = meta[channel_col].values
        markers = None
        if marker_col is not None:
            markers = meta[marker_col].values
        self.config(channels=channels, markers=markers)

        return self

    def config(self, channels=None, markers=None):
        """
        Channel Name: Capitalize abbrivated element name follow with mass number like "Yb137"
        """
        self._channels = channels
        if markers is not None:
            if len(channels) != len(markers):
                print("Unmatched input")
                return self
            markers_map = dict(zip(channels, markers))
            self._markers = OrderedDict((c, markers_map[c]) for c in channels)

        return self

    @property
    def channels(self):
        return self._channels

    @property
    def markers(self):
        return self._markers

    def to_anndata(self):
        # TODO: add multiprocessing support
        lc = len(self._channels)
        lm = len(self._markers)

        if lc == 0:
            try:
                self._channels = read_ROI(self._tree[0]).channels
            finally:
                self._var = pd.DataFrame({'Channels': self._channels})
        elif (lc > 0) & (lm > 0):
            self._var = pd.DataFrame({'Channels': self._channels, 'Markers': list(self._markers.values())})
        elif (lc > 0) & (lm == 0):
            self._var = pd.DataFrame({'Channels': self._channels})
        # anndata require str index, hard set everything to str
        self._var.index = [str(i) for i in range(0, len(self._channels))]
        # print(self._var)

        X = 0
        cells_X = 0
        ann_obs = 0
        for i, d in enumerate(self._tree):
            if len(self._markers) >= 1:
                roi = read_ROI(d).config(channels=self._channels, markers=self._markers)
                exp, cells = roi.exp_matrix()
            else:
                roi = read_ROI(d).config(channels=self._channels)
                exp, cells = roi.exp_matrix()
            cell_count = len(cells)
            obs = np.repeat(np.array([self._obs[i]]), cell_count, axis=0)
            if i == 0:
                X = np.array(exp, dtype=float)
                cells_X = np.array(cells, dtype=object)
                ann_obs = obs
            else:
                X = np.concatenate((X, np.array(exp, dtype=float)), axis=0)
                cells_X = np.concatenate((cells_X, np.array(cells, dtype=object)), axis=0)
                ann_obs = np.concatenate((ann_obs, obs), axis=0)
            print(f"Added: {' '.join(self._obs[i])}")
        # anndata require str index, hard set to str
        ann_obs = pd.DataFrame(ann_obs, columns=self._obs_name, index=[str(i) for i in range(0, len(ann_obs))])
        ann_obs['cell_shape'] = cells_X

        return ad.AnnData(X, obs=ann_obs, var=self._var, dtype=float)


class read_ROI:
    """
    Your .tif/.tiff file should be exported from MCD viewer
    or you can specific channel name in 'page_name' field.
    """

    def __init__(self, folder, mask_pattern="*mask*"):
        """
        Specific the mask image, or automatically select the img name contain "mask".
        """
        self._work_dir = Path(folder)

        # try to find mask
        mask = [i for i in self._work_dir.glob(mask_pattern)]
        try:
            assert len(mask) == 1
        except Exception as e:
            if len(mask) == 0:
                print("No mask found")
                return None
            else:
                print("Found more than one mask")
        self._mask_img = mask[0]

        # get channels info
        self._channels = list()
        self._markers = dict()
        self._channels_files = dict()
        for img in Path(folder).iterdir():
            if img not in mask:
                with tifffile.TiffFile(str(img)) as tif:
                    col_name = tif.pages[0].tags['page_name'].value.decode()
                    pattern = re.compile(r'([a-zA-Z]+)([0-9]{2,})')
                    col = re.findall(pattern, col_name)
                    correct_isotopes = False
                    for c in col:
                        if c[0] in ISOTOPES_NAME:
                            if int(c[1]) in ISOTOPES_MASS_NUMBER_MAP[c[0]]:
                                cname = c[0] + c[1]
                                self._channels.append(cname)
                                self._channels_files[cname] = img
                                correct_isotopes = True
                                break
                    if not correct_isotopes:
                        print(f"Your channel isotope not exists. File: {img}")

    def config_file(self, metadata, channel_col=None, marker_col=None, sep=','):
        meta = pd.read_csv(metadata, sep=sep)
        channels = meta[channel_col].values
        markers = None
        if marker_col is not None:
            markers = meta[marker_col].values
        self.config(channels=channels, markers=markers)

        return self

    def config(self, channels=None, markers=None):
        """
        Channel Name: Capitalize abbrivated element name follow with mass number like "Yb137"
        """
        # TODO: add type check
        selected_channels = []
        not_found_channels = []
        for i, c in enumerate(channels):
            if c in self._channels:
                if c not in selected_channels:
                    selected_channels.append(c)
            else:
                not_found_channels.append(i)
                print(f'{c} not found')
        self._channels = selected_channels
        if markers is not None:
            if len(channels) != len(markers):
                print("Unmatched input")
                return self
            markers_map = dict(zip(channels, markers))
            self._markers = OrderedDict((c, markers_map[c]) for c in selected_channels)

        return self

    @property
    def channels(self):
        return self._channels

    @property
    def markers(self):
        return self._markers

    def exp_matrix(self, polygonize=True):
        if len(self._markers) == 0:
            print("ATTENTION: NO marker specific, using channels' name instead.")
        cells = mask2cells(self._mask_img)

        stacks = np.array([imread(self._channels_files[c]) for c in self._channels])
        data = get_cell_exp_stack(stacks, cells)
        print(f'Detected {data.shape[0]} cells.')
        if polygonize:
            return data, polygonize_cells(cells)
        else:
            return data, cells
