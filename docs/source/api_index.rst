API Reference
=============

There are few modules in spatialtis

- preprocessing
- statistic
- spatial
- plotting

Preprocessing
--------------

.. autoclass:: spatialtis.read_ROIs

    .. automethod:: spatialtis.read_all_ROIs.config
    .. automethod:: spatialtis.read_all_ROIs.config_file
    .. automethod:: spatialtis.read_all_ROIs.to_anndata

Statistic
---------

.. autofunction:: spatialtis.sta.cell_components
.. autofunction:: spatialtis.sta.cell_co_occurrence
.. autofunction:: spatialtis.sta.cell_density
.. autofunction:: spatialtis.sta.cell_morphology


Spatial
-------

.. autoclass:: spatialtis.Neighbors

    .. automethod:: spatialtis.Neighbors.find_neighbors
    .. automethod:: spatialtis.Neighbors.export_neighbors
    .. automethod:: spatialtis.Neighbors.read_neighbors
    .. automethod:: spatialtis.Neighbors.neighbors_count

.. autofunction:: spatialtis.spatial.neighborhood_analysis
.. autofunction:: spatialtis.spatial.spatial_enrichment_analysis
.. autofunction:: spatialtis.spatial.spatial_distribution
.. autofunction:: spatialtis.spatial.spatial_heterogeneity
.. autofunction:: spatialtis.spatial.hotspot
.. autofunction:: spatialtis.spatial.communities


Plotting
--------

Wrapper plotting functions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. autofunction:: spatialtis.plotting.cell_components
.. autofunction:: spatialtis.plotting.cell_co_occurrence
.. autofunction:: spatialtis.plotting.cell_density
.. autofunction:: spatialtis.plotting.cell_morphology
.. autofunction:: spatialtis.plotting.neighborhood_analysis
.. autofunction:: spatialtis.plotting.spatial_enrichment_analysis
.. autofunction:: spatialtis.plotting.spatial_distribution
.. autofunction:: spatialtis.plotting.spatial_heterogeneity


Spatialtis Plotting functions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. autofunction:: spatialtis.plotting.stacked_bar
.. autofunction:: spatialtis.plotting.violin_plot
.. autofunction:: spatialtis.plotting.cell_map
.. autofunction:: spatialtis.plotting.stacked_kde
.. autofunction:: spatialtis.plotting.heatmap
.. autofunction:: spatialtis.plotting.view_palette


Spatialtis color functions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. autofunction:: spatialtis.plotting.get_colors
.. autofunction:: spatialtis.plotting.get_linear_colors


Utility
--------

.. autofunction:: spatialtis.utils.prepare_svca
.. autofunction:: spatialtis.utils.df2adata_uns
.. autofunction:: spatialtis.utils.adata_uns2df
.. autofunction:: spatialtis.utils.col2adata_obs

