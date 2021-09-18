.. automodule:: spatialtis

API Reference
========================

Config
-------
.. autosummary::
    :toctree: api
    :template: autosummary/config_class.rst

    config._Config


IO/Preprocessing
-----------------
.. autosummary::
    :toctree: api
    :nosignatures:

    read_ROIs
    get_result
    transform_points
    transform_shapes


Basic analysis
----------------
.. autosummary::
    :toctree: api
    :nosignatures:

    cell_components
    cell_density
    cell_morphology
    cell_co_occurrence


Spatial analysis
-------------------
.. autosummary::
    :toctree: api
    :nosignatures:

    find_neighbors
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



Base class
-------------

.. autosummary::
    :toctree:
    :nosignatures:

    abc.AnalysisBase


.. currentmodule:: spatialtis.plotting

Plotting
---------

