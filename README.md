<p align="center">
<img src="img/Logo.svg" width="200"/>
<p/>


# SpatialTis
[![Documentation Status](https://readthedocs.org/projects/spatialtis/badge/?version=latest&style=flat-square)](https://spatialtis.readthedocs.io/en/latest/?badge=latest)
![CI](https://flat.badgen.net/github/status/Mr-Milk/SpatialTis?icon=github&label=CI)
![codecov](https://flat.badgen.net/codecov/c/github/Mr-Milk/SpatialTis)
![pypi](https://flat.badgen.net/pypi/v/spatialtis?color=blue)
![licence](https://flat.badgen.net/github/license/Mr-Milk/SpatialTis)

SpatialTis is a high-performance spatial analysis toolkit for single-cell multiplexed tissue data using [`AnnData`](https://icb-anndata.readthedocs-hosted.com/en/stable/#) object as input with **parallel processing** support.

**Documentation**: [![rtd](img/view_on_rtd.svg)](https://spatialtis.readthedocs.io/en/latest/)

[Quick Start](https://spatialtis.readthedocs.io/en/latest/tutorial/1-basic_usage.html)

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
- **Basic analysis**
    - Cell components
    - Cell density
    - Cell morphology
    - Cell co-occurrence
- **Spatial analysis**
    - Find cell neighbors
    - Spatial distribution
    - Spatial heterogeneity
    - Hotspot detection
    - Cell-cell interaction
    - Spatial co-expression
    - Spatial community detection
    - Neighbor dependent markers
