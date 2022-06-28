Getting Started
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

A recommended way to import SpatialTis

    >>> import spatialtis as st
    >>> import spatialtis.plotting as sp

Let's load some data for analysis, SpatialTis takes `AnnData` as input,

    >>> from anndata import read_h5ad
    >>> data = read_h5ad('seqFISH.h5ad')
    >>> data
    AnnData object with n_obs × n_vars = 913 × 9566
        obs: 'Field of View', 'centroid', 'Region', 'cell_type'
        var: 'markers'

We could start to construct the neighbors network and profile cell-cell interactions.

    >>> st.find_neighbors(data,
    ...                   r=180,
    ...                   exp_obs=['Region', 'Field of View'],
    ...                   centroid_key='centroid')
    ⏳ Find neighbors
    🛠 Method: kdtree
    📦 Added to AnnData, obs: 'cell_neighbors'
    📦 Added to AnnData, obs: 'cell_neighbors_count'
    ⏱ 10ms

    >>> cci = st.cell_interactions(data,
    ...                            exp_obs = ['Region', 'Field of View'],
    ...                            centroid_key='centroid',
    ...                            cell_type_key='cell_type')
    ⏳ Cell interaction
    📦 Added to AnnData, uns: 'cell_interaction'
    ⏱ 164ms

Use Config to avoid boilerplate
--------------------------------

Notice that we repeatedly enter :code:`exp_obs` and :code:`centroid_key`. Obviously,
we don't want to write these for every analysis. SpatialTis allows you to set it via Global config.

    >>> from spatialtis import Config
    >>> Config.exp_obs = ['Region', 'Field of View']
    >>> Config.cell_type_key = 'cell_type'
    >>> Config.centroid_key = 'centroid'
    >>> Config.marker_key = 'markers'
    >>> Config.dumps(data)  # save your config
    >>> Config.loads(data)  # load config

Now we could run the analysis in a much cleaner way.

    >>> st.find_neighbors(data)
    >>> st.cell_interaction(data)

Save and get results
---------------------

To access the result of `cell-cell interactions analysis`,
call the :code:`get_result` to get it from `AnnData.uns`

    >>> result = st.get_result(data, 'cell_interaction')

To visualize the analysis.

    >>> sp.cell_interaction(data)
