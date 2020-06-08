Config
--------

SpatialTis allows you to set global configuration, so that you don't have to specify
some parameters in every function.

To set the config::

    from spatialtis import CONFIG

Here is some useful configurations you will deal with a lot, all the names are *CAPITALIZED*.

To let spatialtis how your experiments are designed, this should be a list of `key` from :code:`anndata.obs`.
The order of the name do matters. The last element is assumed to be the ROI level.

    - :code:`CONFIG.EXP_OBS`: :code:`None`

To let spatialtis know where to find your cell type, this should be a `key` name from :code:`anndata.obs`.

    - :code:`CONFIG.CELL_TYPE_COL`: :code:`None`

To set your working environment and OS. Setting working environment to :code:`None` will auto disable progress bar.

    - :code:`CONFIG.WORKING_ENV`: :code:`jupyter` ['jupyter', 'zepplin', None]
    - :code:`CONFIG.OS`: :code:`None` ['Linux', 'Darwin', 'Windows']

To enable paralle processing globally in spatialtis. We used Ray to support this features, which is currently no available
on Windows platform. But you can try WSL2. (This feature won't work if you run spatialtis on Windows)

    - :code:`CONFIG.MULTI_PROCESSING`: :code:`False`

To allocate how many CPU used in paralle processing, the defualt is to use all available resources.

    - :code:`CONFIG.CPU_USED`: :code:`None`

To turn on/off progress bar or you want to config it's appearance.

    - :code:`CONFIG.PROGRESS_BAR`: :code:`True`
    - :code:`CONFIG.PBAR_FORMAT`: :code:`"%s{l_bar}%s{bar}%s{r_bar}%s" % (Fore.GREEN, Fore.CYAN, Fore.GREEN, Fore.RESET,)`


Information stored in the anndata object

    - :code:`CONFIG.CENTROID_COL`: :code:`"centroid"`
    - :code:`CONFIG.COMMUNITY_COL`: :code:`"communities"`
    - :code:`CONFIG.NEIGHBORS_COL`: :code:`"cell_neighbors"`
    - :code:`CONFIG.NEIGHBORS_COUNT_COL`: :code:`"neighbors_count"`
    - :code:`CONFIG.AREA_COL`: :code:`"area"`
    - :code:`CONFIG.SHAPE_COL`: :code:`"cell_shape"`
    - :code:`CONFIG.ECCENTRICITY_COL`: :code:`"eccentricity"`
    - :code:`CONFIG.MARKER_COL`: :code:`"markers"`
    - :code:`CONFIG.CHANNEL_COL`: :code:`"channels"`


