<p align="center">
<img src="https://raw.githubusercontent.com/Mr-Milk/SpatialTis/master/src/favicon-readme.png" width=10%/>
<p/>


# SpatialTis
[![Documentation Status](https://readthedocs.org/projects/spatialtis/badge/?version=latest)](https://spatialtis.readthedocs.io/en/latest/?badge=latest) ![CI](https://github.com/Mr-Milk/SpatialTis/workflows/CI/badge.svg) ![codecov](https://codecov.io/gh/Mr-Milk/SpatialTis/branch/master/graph/badge.svg?token=DYNZ45IPSQ) ![licence](https://badgen.net/badge/licence/Apache%202%2E0/blue)

SpatialTis is a spatial analysis toolkit for single-cell multiplexed tissue data using [`AnnData`](https://icb-anndata.readthedocs-hosted.com/en/stable/#) object as input. Multiprocessing is supported for most analysis functions.

> **Note**: This package is still in development, API may change in the future.

**Documentation**: [![rtd](https://badgen.net/badge/view%20on/read%20the%20docs/blue)](https://spatialtis.readthedocs.io/en/latest/)

**Tutorial**: [![nbviewer](https://badgen.net/badge/view%20on/nbviewer/orange)](https://nbviewer.jupyter.org/github/Mr-Milk/SpatialTis-Tutorial/blob/master/Tutorial%20%28MIBI-dataset%29.ipynb)



## Installation

### Pypi

```shell
pip install spatialtis
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
    - Neighborhood analysis
    - Spatial community detection
    - Marker influence by neighbor cell/marker
