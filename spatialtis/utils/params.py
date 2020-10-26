import functools
import inspect
from typing import Sequence

from spatialtis.config import CONFIG


def get_default_params(func):
    @functools.wraps(func)
    def wrapper(*args, **kw):
        # to handle default parameters
        sig = inspect.signature(func)
        params_names = sig.parameters.keys()
        kw_names = kw.keys()

        if "groupby" in kw_names:
            groupby = kw["groupby"]
            if isinstance(groupby, Sequence):
                kw["groupby"] = list(groupby)
            elif isinstance(groupby, str):
                kw["groupby"] = [groupby]

        def _handle_kwargs(
            kwargs_writer, kwargs, default, error=True, error_message=None
        ):
            if kwargs in params_names:
                if kwargs not in kw_names:
                    if default is not None:
                        kwargs_writer[kwargs] = default
                    else:
                        if error:
                            raise ValueError(error_message)

        # groupby
        _handle_kwargs(
            kw,
            "groupby",
            CONFIG.EXP_OBS,
            True,
            "Experiment observation unclear, set `spatialtis.CONFIG.EXP_OBS` or use argument `groupby=`",
        )
        _handle_kwargs(kw, "mp", CONFIG.MULTI_PROCESSING, False)

        # handle default key
        for a, b, c, d in zip(
            ["type_key", "centroid_key", "metric_key", "shape_key", "marker_key"],
            [
                CONFIG.CELL_TYPE_KEY,
                CONFIG.CENTROID_KEY,
                CONFIG.ECCENTRICITY_KEY,
                CONFIG.SHAPE_KEY,
                CONFIG.MARKER_KEY,
            ],
            ["cell type", "centroid", "metric", "shape", "marker"],
            [
                "CONFIG.CELL_TYPE_KEY",
                "CONFIG.CENTROID_KEY",
                "CONFIG.ECCENTRICITY_KEY",
                "CONFIG.SHAPE_KEY",
                "CONFIG.MARKER_KEY",
            ],
        ):
            _handle_kwargs(
                kw, a, b, True, f"Either specific {c} key using `{a}=` or `{d}`"
            )
        return func(*args, **kw)

    return wrapper
