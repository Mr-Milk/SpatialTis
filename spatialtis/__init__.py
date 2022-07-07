from ._version import version_tuple

__version__ = ".".join([str(i) for i in version_tuple[:3]])

import logging

logging.getLogger("lightning").setLevel(logging.ERROR)
