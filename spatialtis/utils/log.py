import functools
import logging
from time import time

from colorama import Fore

from spatialtis.config import CONFIG

logger = logging.getLogger(__name__)


def pretty_time(t):
    if t < 1:
        return f"{int(t * 1000)}ms"
    elif 1 <= t < 60:
        ms = t * 1000 % 1000
        return f"{int(t)}s{int(ms)}ms"
    elif 60 <= t < 3600:
        minute = t // 60
        second = t % 60
        return f"{int(minute)}m{int(second)}s"
    elif t >= 3600:
        hour = t // 3600
        minute = (t - (hour * 3600)) // 60
        second = t % 60
        return f"{int(hour)}h{int(minute)}m{int(second)}s"


def log_print(text, color="green", verbose=None):
    color_map = {"green": Fore.GREEN, "red": Fore.RED}

    if verbose is None:
        verbose = CONFIG.VERBOSE.INFO

    if CONFIG.WORKING_ENV is None:
        msg = f"{text}"
    else:
        msg = f"{color_map[color]}{text}{Fore.RESET}"
    if verbose:
        logger.info(msg)


def timer(prefix=None, suffix=None, verbose=None):
    """
    Timer decorator to measure the time a function used

    Args:
        prefix: content to add at the front
        suffix: content to add at the tail
        verbose: whether to print

    Returns:

    """
    if verbose is None:
        verbose = CONFIG.VERBOSE.INFO

    def timeit(func):
        @functools.wraps(func)
        def timed(*args, **kw):
            if (prefix is not None) & verbose:
                log_print(prefix, verbose=verbose)
            ts = time()
            result = func(*args, **kw)
            te = time()
            if (suffix is not None) & verbose:
                log_print(suffix, verbose=verbose)
            if verbose:
                if CONFIG.WORKING_ENV is not None:
                    logger.info(
                        f"{Fore.GREEN}Finished! Used {Fore.CYAN}{pretty_time(te - ts)}{Fore.RESET}"
                    )
                else:
                    logger.info(f"Finished! Used {pretty_time(te - ts)}")
            return result

        return timed

    return timeit
