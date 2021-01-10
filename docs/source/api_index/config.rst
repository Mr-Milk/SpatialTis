Config
--------

SpatialTis allows you to set global configuration, so that you don't have to specify
some parameters in every function.

To set the config::

    from spatialtis import CONFIG

Here is some useful configurations you will deal with a lot, all the names are *CAPITALIZED*.


+---------------------------------+-------------------+
| :code:`CONFIG.EXP_OBS`          | :code:`None`      |
+---------------------------------+-------------------+
| :code:`CONFIG.ROI_KEY`          | :code:`None`      |
+---------------------------------+-------------------+
| :code:`CONFIG.CELL_TYPE_KEY`    | :code:`None`      |
+---------------------------------+-------------------+
| :code:`CONFIG.WORKING_ENV`      | :code:`"jupyter"` |
+---------------------------------+-------------------+
| :code:`CONFIG.OS`               | :code:`None`      |
+---------------------------------+-------------------+
| :code:`CONFIG.VERBOSE`          | :code:`True`      |
+---------------------------------+-------------------+
| :code:`CONFIG.MULTI_PROCESSING` | :code:`False`     |
+---------------------------------+-------------------+

Storage keys, please explicitly specify for your own data

+---------------------------------+------------------------+
| :code:`CONFIG.CENTROID_KEY`     | :code:`"centroid"`     |
+---------------------------------+------------------------+
| :code:`CONFIG.SHAPE_KEY`        | :code:`"cell_shape"`   |
+---------------------------------+------------------------+
| :code:`CONFIG.AREA_KEY`         | :code:`"area"`         |
+---------------------------------+------------------------+
| :code:`CONFIG.ECCENTRICITY_KEY` | :code:`"eccentricity"` |
+---------------------------------+------------------------+
| :code:`CONFIG.MARKER_KEY`       | :code:`"markers"`      |
+---------------------------------+------------------------+


CONFIG.EXP_OBS
=================

To let spatialtis how your experiments are designed, this should be a list of `keys` from :code:`anndata.obs`.
The order of the name do matters. The last element is assumed to be the ROI level. If not, you must specify using
:code:`CONFIG.ROI_KEY`.


CONFIG.CELL_TYPE_KEY
=====================
To let spatialtis know where to find your cell type, this should be a `key` name from :code:`anndata.obs`.


CONFIG.WORKING_ENV
===================

Available options: :code:`["jupyter", "zepplin", None]`.
To set working environment to :code:`None` will automatically disable progress bar.


CONFIG.OS
===================

Available options: :code:`["Linux", "Darwin", "Windows"]`.
Normally, you don't need to specific your system, it will be auto detected by SpatialTis


CONFIG.MULTI_PROCESSING
=========================

To enable paralle processing globally in spatialtis. We used Ray to support this features, the support for windows platform
is still in experimental stage, if something went wrong, you should turn it off.


CONFIG.VERBOSE
===============

Config the amount of information print by SpatialTis

:code:`VERBOSE.ANNDATA`: :code:`[True, False]` whether to display the change made to anndata object

:code:`VERBOSE.PBAR`: :code:`[True, False]` whether to display the progress bar.

:code:`VERBOSE.INFO`: :code:`[True, False]` whether to display the runtime info and timer.

CONFIG.PBAR_FORMAT
====================

Default value: :code:`"%s{l_bar}%s{bar}%s{r_bar}%s" % (Fore.GREEN, Fore.CYAN, Fore.GREEN, Fore.RESET,)`

To configure the appearance of progress bar

