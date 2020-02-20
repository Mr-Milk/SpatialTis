import platform

from .preprocessing import read_all_ROIs, read_ROI

if platform.system() in ["Linux", "Darwin"]:
    import logging
    import ray

    ray.init(logging_level=logging.FATAL, ignore_reinit_error=True)