import re
from collections import OrderedDict
from pathlib import Path
from typing import Sequence, Union

import numpy as np
import pandas as pd
from skimage.external import tifffile
from skimage.io import imread

from spatialtis.config import ISOTOPES_MASS_NUMBER_MAP, ISOTOPES_NAME

from ._geom import geom_cells
from ..config import CONFIG


def mask2cells(mask_img: Union[Path, str], ignore_bg: bool = True) -> Sequence:
    """
    Parameters
    mask_img: Path/str, the path to mask img
    ignore_bg: bool, ignore background in mask
    """
    # read mask image
    mask = imread(mask_img)
    # number of cells
    counts = np.unique(mask)
    # in case the number in mask is not continuous
    mapper = dict(zip(counts, range(0, len(counts))))
    # create list for each cell
    cells = [[] for i in range(0, len(counts))]
    # find the exact points that belong to each cell
    iter = np.nditer(mask, flags=["multi_index"])
    while not iter.finished:
        cells[mapper[int(iter[0])]].append(iter.multi_index)
        iter.iternext()
    # usually 0 for background
    # and if index 0 of cells didn't contain only one cell
    if ignore_bg:
        return cells[1:]
    else:
        return cells


def get_cell_exp_stack(
        stack: Sequence, cells: Sequence, method: str = "mean"
) -> Sequence:
    """
    Parameters
    channel: list or ndarray, matrix info of the channel
    cells: list or ndarray, each element contains points for each cell
    method: str, ('mean' / 'median' / 'sum' / ...)  any numpy method to compute expression level of single cell
    """
    cells_density = list()
    for cell in cells:
        cell_pixels = stack[:, [i[0] for i in cell], [i[1] for i in cell]]
        exec(f"cells_density.append([np.{method}(density) for density in cell_pixels])")
    # each secondary array as a channel
    return cells_density


class ROIreader:
    tree = None
    markers = None
    channels = None
    var = None
    _var = None

    def config_(self, channels=None, markers=None, callback=None):
        """
        Channel Name: Capitalize abbreviated element name follow with mass number like "Yb137"
        """
        if markers is not None:
            if len(channels) != len(markers):
                # TODO: fix uninformative print
                print("Unmatched input")
                return self
            markers_map = dict(zip(channels, markers))
            self.markers = OrderedDict((c, markers_map[c]) for c in channels)

        if callback is not None:
            try:
                callback(self)
            except NameError:
                print("callback is not a function")

        return self

    def config_file_(
            self, metadata, channel_col=None, marker_col=None, sep=",", callback=None
    ):
        meta = pd.read_csv(metadata, sep=sep)
        channels = meta[channel_col].values
        markers = None
        if marker_col is not None:
            markers = meta[marker_col].values
        self.config_(channels=channels, markers=markers)

        if callback is not None:
            try:
                callback(self)
            except NameError:
                print("callback is not a function")

        return self

    def set_info(self):
        lc = len(self.channels)
        lm = len(self.markers)
        channel_key = CONFIG.CHANNEL_COL
        marker_key = CONFIG.MARKER_COL

        if lc == 0:
            try:
                self.channels = read_ROI(self.tree[0]).channels
            finally:
                self._var = pd.DataFrame({channel_key: self.channels})
        elif (lc > 0) & (lm > 0):
            self.var = pd.DataFrame(
                {channel_key: self.channels, marker_key: list(self.markers.values())}
            )
        elif (lc > 0) & (lm == 0):
            self.var = pd.DataFrame({marker_key: self.channels})
        # anndata require str index, hard set everything to str
        self.var.index = [str(i) for i in range(0, len(self.channels))]


class read_ROI(ROIreader):
    """
    Your .tif/.tiff file should be exported from MCD viewer
    or you can specific channel name in 'page_name' field.
    """

    def __init__(self, folder, mask_pattern="*mask*", stacked=False):
        """
        Specific the mask image, or automatically select the img name contain "mask".
        """
        self.__stacked = stacked
        self.__work_dir = Path(folder)

        # try to find mask
        mask = [i for i in self.__work_dir.glob(mask_pattern)]
        if len(mask) == 0:
            raise FileNotFoundError(f"Mask not found in {str(folder)}")
        elif len(mask) > 1:
            raise ValueError(f"Found more than one mask in {str(folder)}")
        self.__mask_img = mask[0]

        # get channels info
        self.channels = list()
        self.markers = dict()
        self.__channels_files = dict()
        self.__stacks = 0

        if stacked:
            stacks_count = 0
            for img in Path(folder).iterdir():
                if img not in mask:
                    # From skimage doc: The different color bands/channels are stored in the third dimension
                    # so we need to transpose it
                    self.__stacks = np.transpose(imread(str(img)), (2, 1, 0))
                    stacks_count += 1

            if stacks_count == 0:
                raise FileNotFoundError(f"ROI image not found in {str(folder)}")
            elif stacks_count > 1:
                raise ValueError(f"Found more than one ROI image in {str(folder)}")

        else:
            for img in Path(folder).iterdir():
                if img not in mask:
                    with tifffile.TiffFile(str(img)) as tif:
                        col_name = tif.pages[0].tags["page_name"].value.decode()
                        pattern = re.compile(r"([a-zA-Z]+)([0-9]{2,})")
                        col = re.findall(pattern, col_name)
                        correct_isotopes = False
                        for c in col:
                            if c[0] in ISOTOPES_NAME:
                                if int(c[1]) in ISOTOPES_MASS_NUMBER_MAP[c[0]]:
                                    cname = c[0] + c[1]
                                    self.channels.append(cname)
                                    self.__channels_files[cname] = img
                                    correct_isotopes = True
                                    break
                        if not correct_isotopes:
                            print(f"Your channel isotope not exists. File: {str(img)}")

    def config(self, channels=None, markers=None):

        # selected_channels = filter_channels(self, channels=channels)
        super().config_(channels=channels, markers=markers)
        if not self.__stacked:
            self.__stacks = np.asarray(
                [
                    tifffile.TiffFile(str(self.__channels_files[c])).pages[0].asarray()
                    for c in self.channels
                ]
            )
        return self

    def exp_matrix(self, method="mean", polygonize="convex", alpha=0):
        if len(self.markers) == 0:
            print("ATTENTION: NO marker specific, using channels' name instead.")
        cells = mask2cells(self.__mask_img)
        cells, geom_info = geom_cells(cells, method=polygonize, alpha=alpha)
        data = get_cell_exp_stack(self.__stacks, cells, method=method)
        # print(f"Detected {len(data)} cells.")

        return data, geom_info
