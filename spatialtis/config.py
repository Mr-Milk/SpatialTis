"""
Setting Global config for whole processing level
"""
from typing import Optional, Sequence, Union

# can be override in plotting function
# ['jupyter', 'zepplin', None]
WORKING_ENV: Optional[str] = 'jupyter'

"""
If working in Jupyter Lab/Hub, please install the following dependencies
It should be included in the virtual environment for spatialTis

    ```jupyter labextension install @bokeh/jupyter_bokeh```

If working with Zeppelin, please install `bkzep` (Only PyPI)

    ```pip install bkzep```

"""

# Which file format to save across spatialTis
# ['png', 'svg', 'html', None]
SAVE_FORMAT: Union[str, Sequence[str], None] = None
"""
Supported format:
Static file: PNG, SVG,
Interactive file: HTML
The default will embed the JS resource inside html
if you want to switch to CDN, set `use_CDN = True`
"""

# Embed bokeh JS resource or use CDN
use_CDN: bool = False

ROI_KEY: str = 'ROI'
