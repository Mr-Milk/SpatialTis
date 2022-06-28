Input Data
===============

SpatialTis use :code:`AnnData` as input format.

ROI Information
---------------

If you want to process multiple ROIs at the same time,
It's important that you annotate your cell with ROI.

An example table from :code:`AnnData.obs` will look like this:

.. csv-table::
    :file: ../public/roi_table.csv
    :header-rows: 1

In this table, the columns **Stage**, **Patient**, **Tissue** and **ROI** tells
the information of ROI. Those information wil be shown in the visualization.
The column `ROI` should assign a unique name for each ROI.


Spatial Coordination
--------------------

Centroid
++++++++++

SpatialTis accepts different kinds of spatial coordination format.

- (Default) One key in :code:`.obsm`, it should be a :code:`numpy.ndarray`, default will read from **spatial**
- Two/three columns in :code:`.obs`, each store one dimension of coordination.
  Respectively, `x`, `y` and `z` if you have 3D data.
- One column in :code:`.obs` stored as
  `wkt <https://en.wikipedia.org/wiki/Well-known_text_representation_of_geometry>`_
  (well-known text) format, which could be easily serialized and deserialized.
  Or python list/tuple.

If you have two columns `X` and `Y` in `AnnData.obs` that store x, y coordination,
you could transform them into wkt format:

>>> import spatialtis as st
>>> st.wkt_points(data, ('X', 'Y'))

Or if you save your coordination in single column :code:`centroid`, transform it like:

>>> st.wkt_points(data, 'centroid')

Cell Shape (Polygon)
+++++++++++++++++++++

SpatialTis only reads polygons in wkt format, if you don't have wkt ready format.
You should store a polygon in an array container like
:code:`[(1, 2), (3, 4), ..., (100, 100)]`. Store all your polygons in one column.
You can now transform it to wkt format:

>>> st.wkt_shapes(data, 'shape')

