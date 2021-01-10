.. spatialTis-doc documentation master file, created by
   sphinx-quickstart on Wed Jan 22 08:48:16 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to SpatialTis's documentation!
==========================================

.. image:: img/Logo.svg
   :align: center
   :width: 20%

|doc| |ci| |coverage| |pypi| |license|

.. |doc| image:: https://readthedocs.org/projects/spatialtis/badge/?version=latest&style=flat-square
.. |ci| image:: https://flat.badgen.net/github/status/Mr-Milk/SpatialTis?icon=github&label=CI
.. |coverage| image:: https://flat.badgen.net/codecov/c/github/Mr-Milk/SpatialTis
.. |license| image:: https://flat.badgen.net/github/license/Mr-Milk/SpatialTis
.. |pypi| image:: https://flat.badgen.net/pypi/v/spatialtis?color=blue

SpatialTis is a spatial analysis toolkit for single-cell multiplexed tissue data.
For a quick view of what it's does, check our following examples:

- MIBI breast cancer 180K cells: |MIBI-tutorial|_
- IMC diabetes 1.7M cells: |IMC-tutorial|_

.. |MIBI-tutorial| image:: ../../img/view_on_nbviewer.svg
.. _MIBI-tutorial: https://nbviewer.jupyter.org/github/Mr-Milk/SpatialTis-Tutorial/blob/master/Tutorial-1%20%28MIBI-dataset%29.ipynb

.. |IMC-tutorial| image:: ../../img/view_on_nbviewer.svg
.. _IMC-tutorial: https://nbviewer.jupyter.org/github/Mr-Milk/SpatialTis-Tutorial/blob/master/Tutorial-2%20%28IMC-dataset%29.ipynb


.. toctree::
   :maxdepth: 1
   :glob:
   :caption: Usage

   usage/get_started
   usage/example
   usage/installation
   usage/image_export
   usage/notebook

.. toctree::
   :maxdepth: 1
   :glob:
   :caption: Tutorial

   tutorial/basic_usage
   tutorial/set_config
   tutorial/plotting
   tutorial/example1
   tutorial/example2

.. toctree::
   :maxdepth: 1
   :glob:
   :caption: About

   api_index/api_index
   about/implementation
   about/development_guide
   about/QA


