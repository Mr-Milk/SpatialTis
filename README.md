# SpatialTis
spatialtis is a toolkit for analysis of multiplexed tissue imaging data written in python. Applied to any multiplexed tissue data with single cell solution. eg. Imaging Mass Cytometry and Multiplex Ion Beam Imaging.

spatialtis uses [`annData`](https://icb-anndata.readthedocs-hosted.com/en/stable/#) object as input, to integrated with popular single cell analysis library [`Scanpy`](https://scanpy.readthedocs.io/en/latest/index.html) . Preprocessing module of spatialtis allow you to transfrom your data into `annData` objects, please prepare your multiplexed image data and cell mask for each ROI in seperated folders.
