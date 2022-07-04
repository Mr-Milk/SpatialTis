API Reference
========================

.. module:: spatialtis.config

.. currentmodule:: spatialtis

Config
-------
.. autosummary::
    :toctree: API
    :template: autosummary/config_class.rst

    config._Config


.. module:: spatialtis.preprocessing

.. currentmodule:: spatialtis

Cell quantification from Images
--------------------------------
.. autosummary::
    :toctree: API
    :nosignatures:

    read_images


Read from 10x visium
--------------------------------
.. autosummary::
    :toctree: API
    :nosignatures:

    read_visium


WKT Format helper
--------------------------------
.. autosummary::
    :toctree: API
    :nosignatures:

    wkt_points
    wkt_shapes

IO
--------------------------------
.. autosummary::
    :toctree: API
    :nosignatures:

    get_result

.. module:: spatialtis.basic

.. currentmodule:: spatialtis


Basic analysis
----------------
.. autosummary::
    :toctree: API
    :nosignatures:

    cell_components
    cell_density
    cell_morphology
    cell_co_occurrence


.. module:: spatialtis.spatial

.. currentmodule:: spatialtis

Spatial analysis
-------------------
.. autosummary::
    :toctree: API
    :nosignatures:

    find_neighbors
    spatial_weights
    hotspot
    cell_dispersion
    cell_interaction
    cell_community
    spatial_heterogeneity
    spatial_autocorr
    spatial_enrichment
    spatial_coexp
    NCD_marker
    NMD_marker
    somde
    GCNG

.. module:: spatialtis.abc

.. currentmodule:: spatialtis

Base class
-------------

.. autosummary::
   :toctree: API
   :template: autosummary/abc.rst

    abc.AnalysisBase


.. module:: spatialtis.plotting

.. currentmodule:: spatialtis

To use plotting, import it in following schema:

>>> import spatialtis.plotting as sp


ROI Visualization
------------------

.. autosummary::
    :toctree: API
    :nosignatures:

    plotting.cell_map
    plotting.expression_map


Analysis Visualization
-----------------------

.. autosummary::
    :toctree: API
    :nosignatures:

    plotting.cell_components
    plotting.cell_density
    plotting.cell_morphology
    plotting.cell_co_occurrence
    plotting.cell_dispersion
    plotting.spatial_heterogeneity
    plotting.cell_interaction
    plotting.spatial_enrichment

