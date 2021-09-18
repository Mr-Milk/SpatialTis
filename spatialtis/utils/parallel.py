import logging
from typing import Callable

from .iters import pbar_iter


def create_remote(funcs: Callable):
    try:
        import ray
        if not ray.is_initialized():
            ray.init(
                logging_level=logging.FATAL,
                ignore_reinit_error=True,
            )
    except Exception:
        raise SystemError("Initiate ray failed, parallel processing not available.")

    if isinstance(funcs, Callable):
        return ray.remote(funcs)
    else:
        raise TypeError("Must be a function or a list of function")


def run_ray(jobs, desc=""):
    try:
        import ray
    except Exception:
        raise SystemError("Initiate ray failed, parallel processing not available.")

    def exec_iterator(obj_ids):
        while obj_ids:
            done, obj_ids = ray.wait(obj_ids)
            yield ray.get(done[0])

    for _ in pbar_iter(exec_iterator(jobs), desc=desc, total=len(jobs)):
        pass

    return ray.get(jobs)
