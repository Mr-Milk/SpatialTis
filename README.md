<p align="center">
<img src="https://raw.githubusercontent.com/Mr-Milk/SpatialTis/master/img/Logo.svg" width="200"/>
<p/>


# SpatialTis
[![Documentation Status](https://readthedocs.org/projects/spatialtis/badge/?version=latest&style=flat-square)](https://spatialtis.readthedocs.io/en/latest/?badge=latest)
![CI](https://flat.badgen.net/github/status/Mr-Milk/SpatialTis?icon=github&label=CI)
![codecov](https://flat.badgen.net/codecov/c/github/Mr-Milk/SpatialTis)
![pypi](https://flat.badgen.net/pypi/v/spatialtis?color=blue)
![licence](https://flat.badgen.net/github/license/Mr-Milk/SpatialTis)

SpatialTis is an ultra-fast spatial analysis toolkit for large-scale spatial single-cell data.

- ✔️ Spatial Transcriptome (Non single-cell)
- ✔️ Spatial Proteome (Single-cell)
- 🦀 Core algorithms implements in Rust
- 🚀 Parallel processing support

### 🔋 Highlighted spatial analysis

- Cell neighbors search (KD-Tree/R-Tree/Delaunay)
- Cell-Cell Interaction
- Marker spatial co-expression
- Spatial variable genes (current support: SOMDE)
- GCNG: Inferring ligand-receptor using graph convolution network
- Identify neighbor dependent markers

### 📦 Other analysis

  - Spatial distribution
  - Hotspot detection
  - Spatial auto-correlation
  - Spatial heterogeneity

[Quick Start](https://spatialtis.readthedocs.io/en/latest/usage/get_started.html)


## Installation

### pypi

SpatialTis requires **Python >= 3.8**.

![Version Support](https://img.shields.io/pypi/pyversions/spatialtis?style=flat-square)

```shell
pip install spatialtis

# For full features
pip install 'spatialtis[all]'
```

Install the current development version

```shell
pip install git+https://github.com/Mr-Milk/SpatialTis.git
```

### Docker

```shell
docker pull mrmilk/spatialtis
```

To start a docker container:
```shell
cd your/data/
docker run -it --rm -p 8888:8888 -v "${PWD}:/analysis" spatialtis
```
- `-it`: Run the container in interactive mode
- `-rm`: Clean file system in container after shutting down
- If local port 8888 is taken, try `-p 9999:8888` and change to 9999. 
- `-v`: Mount your data directory to the working directory `/analysis` 
  in the container. `${PWD}` is the directory where you run this command. 
  All changes made in this directory will be saved.


## Low level API

If you are interested in using low level algorithms yourself,
Please refer to [spatialtis_core](https://github.com/Mr-Milk/SpatialTis-core)
It provides clear document for all exposed API.
