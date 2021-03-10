Getting Started!
=================

To install SpatialTis, please see the `Installation <installation.rst>`_ section. It's recommended to start a new virtual environment before installing the SpatialTis.

Prerequisites
-------------

Must learn:

- `anndata <https://anndata.readthedocs.io/en/stable/>`_: `Tutorial <https://falexwolf.de/blog/171223_AnnData_indexing_views_HDF5-backing/>`_ written by the author of `anndata` and `scanpy`

Optional:

- `numpy <https://numpy.org/devdocs/user/quickstart.html>`_: For array operation
- `pandas <https://pandas.pydata.org/pandas-docs/stable/getting_started/tutorials.html>`_: For table operation
- `scanpy <https://scanpy-tutorials.readthedocs.io/en/latest/index.html>`_: Used for cell type calling and scRNA-seq related analysis

Basic Usage
--------------------------

Before any analysis using SpatialTis, it's necessary to set up some :doc:`global config <../api_index/config>`,
so that you don't need to specific them every time::

    from spatialtis import CONFIG

    CONFIG.EXP_OBS = ['Patient', 'Sample', 'ROI']
    CONFIG.ROI_KEY = 'ROI'
    CONFIG.CELL_TYPE_KEY = 'cell_type'
    CONFIG.CENTROID_KEY = 'centroid'
    CONFIG.MARKERS_KEY = 'markers'

There are an analysis module and a visualization module `plotting` in SpatialTis::

    import spatialtis as st
    import spatialtis.plotting as sp

First load the data::

    from anndata import read_h5ad
    data = read_h5ad('sample.h5ad')

Normally an analysis function comes with a visualization function,
they share the same name but registered in different modules.
**DON'T** import these functions individually, it will cause conflicts

Call an analysis function for cell_components::

    st.cell_components(data)

Call visualization function for cell_components::

    sp.cell_components(data)

Now that you've learn some basic of SpatialTis, feel free to play around. If you want to know more about analysis and
visualization with SpatialTis, check the :doc:`Tutorial <tutorial>`.



