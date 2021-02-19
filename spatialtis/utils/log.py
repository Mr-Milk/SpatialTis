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
