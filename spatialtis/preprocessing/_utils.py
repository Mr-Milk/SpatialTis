import numpy as np
import pandas as pd
from pathlib import Path
from skimage.io import imread
from collections import OrderedDict

from typing import Union, Sequence


def mask2cells(
        mask_img: Union[Path, str],
        ignore_bg: bool = True
) -> Sequence:
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
        channel: Sequence,
        cells: Sequence,
        method: str = 'mean'
) -> Sequence:
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
        stack: Sequence,
        cells: Sequence,
        method: str = 'mean'
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
        exec(f'cells_density.append([np.{method}(density) for density in cell_pixels])')
    # each secondary array as a channel
    return cells_density


def config(cls, channels=None, markers=None, callback=None):
    """
    Channel Name: Capitalize abbrivated element name follow with mass number like "Yb137"
    """
    cls._channels = channels
    if markers is not None:
        if len(channels) != len(markers):
            print("Unmatched input")
            return cls
        markers_map = dict(zip(channels, markers))
        cls._markers = OrderedDict((c, markers_map[c]) for c in channels)

    if callback is not None:
        try:
            callback(cls)
        except NameError:
            print("callback is not a function")

    return cls


def config_file(cls, metadata, channel_col=None, marker_col=None, sep=',', callback=None):
    meta = pd.read_csv(metadata, sep=sep)
    channels = meta[channel_col].values
    markers = None
    if marker_col is not None:
        markers = meta[marker_col].values
    cls.config(channels=channels, markers=markers)

    if callback is not None:
        try:
            callback(cls)
        except NameError:
            print("callback is not a function")

    return cls


def filter_channels(cls, channels=None):
    # TODO: add type check
    selected_channels = []
    not_found_channels = []
    for i, c in enumerate(channels):
        if c in cls._channels:
            if c not in selected_channels:
                selected_channels.append(c)
        else:
            not_found_channels.append(i)
            print(f"{c} not found")
    return selected_channels

