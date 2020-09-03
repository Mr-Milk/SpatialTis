.. spatialTis-doc documentation master file, created by
   sphinx-quickstart on Wed Jan 22 08:48:16 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to SpatialTis's documentation!
==========================================

.. image:: src/favicon-chrome.png
   :align: center
   :width: 20%

SpatialTis is a spatial analysis toolkit for single-cell multiplexed tissue data.

For a quick view of what it's does, check our following examples:

- MIBI breast cancer 180K cells: |MIBI-tutorial|_
- IMC diabetes 1.7M cells: |IMC-tutorial|

.. |MIBI-tutorial| image:: https://badgen.net/badge/view%20on/nbviewer/orange
.. _MIBI-tutorial: https://nbviewer.jupyter.org/github/Mr-Milk/SpatialTis-Tutorial/blob/master/Tutorial%20%28MIBI-dataset%29.ipynb

.. |IMC-tutorial| image:: https://badgen.net/badge/view%20on/nbviewer/orange


.. toctree::
   :maxdepth: 2
   :glob:
   :caption: Usage

   usage/get_started
   usage/example
   usage/installation
   usage/image_export

.. toctree::
   :maxdepth: 2
   :glob:
   :caption: About

   api_index/api_index
   about/implementation
   about/development_guide
   about/QA


