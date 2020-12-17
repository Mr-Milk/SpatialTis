import functools
from time import time

from spatialtis.config import CONFIG
from spatialtis.console import console


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


def log_print(text, custom=False, color="green", verbose=None):
    if verbose is None:
        verbose = CONFIG.VERBOSE

    if verbose:
        if not custom:
            console.print(text, style=color)
        else:
            console.print(text)


def timer(task_name=None, suffix=None):
    """
    Timer decorator to measure the time a function used

    Args:
        task_name: content to add at the front
        suffix: content to add at the tail

    Returns:

    """

    def timeit(func):
        @functools.wraps(func)
        def timed(*args, **kw):
            if task_name is not None:
                log_print(
                    f":hourglass_not_done: [green]{task_name}[/green]", custom=True
                )
            ts = time()
            result = func(*args, **kw)
            te = time()
            if suffix is not None:
                log_print(suffix)
            log_print(
                f":stopwatch: [green]Finished![/green] [bold cyan]{pretty_time(te - ts)}[/bold cyan]",
                custom=True,
            )
            return result

        return timed

    return timeit
