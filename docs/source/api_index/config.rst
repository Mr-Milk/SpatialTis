Config
--------

SpatialTis allows you to set global configuration, so that you don't have to specify
some parameters in every function.

To set the config::

    from spatialtis import CONFIG

Here is some useful configurations you will deal with a lot, all the names are *CAPITALIZED*.


+---------------------+---------------------------------------+----------------------+------------------+
| Options             | Description                           | Default Value        | Type             |
+---------------------+---------------------------------------+----------------------+------------------+
| :code:`OS`          | Operating system                      | Detect automatically | `str`            |
+---------------------+---------------------------------------+----------------------+------------------+
| :code:`MP`          | Use multi-processing                  | :code:`True`         | `bool`           |
+---------------------+---------------------------------------+----------------------+------------------+
| :code:`WORKING_ENV` | The current working environment       | Detect automatically | `str`            |
+---------------------+---------------------------------------+----------------------+------------------+
| :code:`VERBOSE`     | Control the print statement           | :code:`True`         | `bool`           |
+---------------------+---------------------------------------+----------------------+------------------+
| :code:`PBAR`        | Control the process bar               | :code:`True`         | `bool`           |
+---------------------+---------------------------------------+----------------------+------------------+
| :code:`AUTO_SAVE`   | To autosave the visualization results | :code:`False`        | `bool` or `Path` |
+---------------------+---------------------------------------+----------------------+------------------+
| :code:`EXP_OBS`     | The design of the experiment data     | :code:`None`         | `List[str]`      |
+---------------------+---------------------------------------+----------------------+------------------+
| :code:`ROI_KEY`     | The key of ROI                        | :code:`None`         | `str`            |
+---------------------+---------------------------------------+----------------------+------------------+

Storage keys, please explicitly specify for your own data

+--------------------------+-----------------------------------+------------------------+-------+
| Options                  | Description                       | Default Value          | Type  |
+--------------------------+-----------------------------------+------------------------+-------+
| :code:`CELL_TYPE_KEY`    | Cell type key in `AnnData.obs`    | :code:`None`           | `str` |
+--------------------------+-----------------------------------+------------------------+-------+
| :code:`MARKER_KEY`       | Marker key in `AnnData.var`       | :code:`marker`         | `str` |
+--------------------------+-----------------------------------+------------------------+-------+
| :code:`NEIGHBORS_KEY`    | Neighbors key in `AnnData.obs`    | :code:`cell_neighbors` | `str` |
+--------------------------+-----------------------------------+------------------------+-------+
| :code:`CENTROID_KEY`     | Centroid key in `AnnData.obs`     | :code:`centroid`       | `str` |
+--------------------------+-----------------------------------+------------------------+-------+
| :code:`AREA_KEY`         | Area key in `AnnData.obs`         | :code:`area`           | `str` |
+--------------------------+-----------------------------------+------------------------+-------+
| :code:`SHAPE_KEY`        | Shape key in `AnnData.obs`        | :code:`cell_shape`     | `str` |
+--------------------------+-----------------------------------+------------------------+-------+
| :code:`ECCENTRICITY_KEY` | Eccentricity key in `AnnData.obs` | :code:`eccentricity`   | `str` |
+--------------------------+-----------------------------------+------------------------+-------+


