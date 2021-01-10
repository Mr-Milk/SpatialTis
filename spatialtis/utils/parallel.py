import logging
from typing import Callable, List, Union

from tqdm import tqdm


def create_remote(funcs: Union[List[Callable], Callable]):
    try:
        import ray

        ray.init(
            logging_level=logging.FATAL, ignore_reinit_error=True,
        )
    except Exception:
        raise SystemError("Initiate ray failed, parallel processing not available.")

    if isinstance(funcs, Callable):
        return ray.remote(funcs)
    elif isinstance(funcs, List):
        return [ray.remote(f) for f in funcs]
    else:
        raise TypeError("Must be a function or a list of function")


def run_ray(jobs, pbar_config=None):
    try:
        import ray
    except Exception:
        raise SystemError("Initiate ray failed, parallel processing not available.")

    def exec_iterator(obj_ids):
        while obj_ids:
            done, obj_ids = ray.wait(obj_ids)
            yield ray.get(done[0])

    for _ in tqdm(exec_iterator(jobs), **pbar_config):
        pass

    results = ray.get(jobs)

    return results
