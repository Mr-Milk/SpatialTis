<p align="center">
<img src="https://raw.githubusercontent.com/Mr-Milk/SpatialTis/master/src/Logo.svg" width="200"/>
<p/>


# SpatialTis
[![Documentation Status](https://readthedocs.org/projects/spatialtis/badge/?version=latest&style=flat-square)](https://spatialtis.readthedocs.io/en/latest/?badge=latest)
![CI](https://flat.badgen.net/github/status/Mr-Milk/SpatialTis?icon=github&label=CI)
![codecov](https://flat.badgen.net/codecov/c/github/Mr-Milk/SpatialTis)
![pypi](https://flat.badgen.net/pypi/v/spatialtis?color=blue)
![licence](https://flat.badgen.net/github/license/Mr-Milk/SpatialTis)

SpatialTis is a high-performance spatial analysis toolkit for single-cell multiplexed tissue data using [`AnnData`](https://icb-anndata.readthedocs-hosted.com/en/stable/#) object as input with **parallel processing** support.

**Documentation**: [![rtd](img/view_on_rtd.svg)](https://spatialtis.readthedocs.io/en/latest/)

**Tutorial**: 
- MIBI Data (Breast cancer, 180K cells) [![nbviewer](img/view_on_nbviewer.svg)](https://nbviewer.jupyter.org/github/Mr-Milk/SpatialTis-Tutorial/blob/master/Tutorial-1%20%28MIBI-dataset%29.ipynb) | [Download data](https://uofmacau-my.sharepoint.com/:u:/g/personal/yb97643_umac_mo/ET7-chqWIc9EqSEtY-foQ7IBURusGw9hlTSBC3xD_bNdgw?download=1)
- IMC Data (Diabetes, 1.7M cells) [![nbviewer](img/view_on_nbviewer.svg)](https://nbviewer.jupyter.org/github/Mr-Milk/SpatialTis-Tutorial/blob/master/Tutorial-2%20%28IMC-dataset%29.ipynb) | [Download data](https://uofmacau-my.sharepoint.com/:u:/g/personal/yb97643_umac_mo/EXJFp1Nn_k5NphOp986lGvABmDNC_fNPGjrw5xN4NUPnRA?download=1)

[Download](https://github.com/Mr-Milk/SpatialTis-Tutorial) the examples and try it on your own, see how fast SpatialTis is.

## Installation

### pypi

Install the basics

```shell
pip install spatialtis
```

For the full features

```shell
pip install 'spatialtis[all]'
```

Install the current development version

```shell
pip install git+https://github.com/Mr-Milk/SpatialTis.git
```

## SpatialTis modules

- **Preprocessing**
- **Data statistic**
    - Cell components
    - Cell density
    - Cell morphology
    - Cell co-occurrence
- **Find cell neighbors**
- **Spatial analysis**
    - Spatial distribution
    - Spatial heterogeneity
    - Hotspot detection
    - Cell-cell interaction
    - Markers co-expression
    - Spatial community detection
    - Neighbor dependent markers
