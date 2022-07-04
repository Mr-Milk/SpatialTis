from rich.progress import track
from typing import Iterable

from spatialtis.config import Config, console


def pbar_iter(obj: Iterable, desc: str = None, **kwargs):
    for i in track(
            obj, disable=(not Config.verbose), console=console, description=desc, **kwargs
    ):
        yield i
