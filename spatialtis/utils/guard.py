from types import ModuleType
from typing import List
from importlib import import_module

from thefuzz import process

from spatialtis.abc import AnalysisBase


class NeighborsNotFoundError(Exception):
    pass


class CellTypeNotFoundError(Exception):
    pass


def check_neighbors(ab: AnalysisBase):
    if not ab.neighbors_exists:
        raise NeighborsNotFoundError("Neighbors not found! Run `spatialtis.find_neighbors` first.")


def check_cell_type(ab: AnalysisBase):
    if not ab.has_cell_type:
        raise CellTypeNotFoundError("Cell Type not found! Please set `cell_type_key`")


def options_guard(select: str, options: List[str]):
    if select not in options:
        possible_select = process.extractOne(select, options)[0]
        raise ValueError(f"'{select}' is not available, do you mean '{possible_select}'. "
                         f"Available options are {', '.join([str(i) for i in options])}")
    else:
        return select


def try_import(module_name: str, install_name=None) -> ModuleType:
    if install_name is None:
        install_name = module_name
    try:
        module = import_module(module_name)
    except (ImportError, ModuleNotFoundError):
        raise ModuleNotFoundError(f"Cannot import {module_name}, try `pip install {install_name}` to install it.")
    return module
