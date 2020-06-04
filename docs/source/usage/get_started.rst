Getting Started!
=================

To install spatialtis, please see the `Installation <installation.rst>`_ section. It's recommended to start a new virtual environment before installing the spatialtis.

Prerequisites
-------------

To work with spatialtis, you need to know about `anndata`, if you are familiar with `pandas` and `numpy` this should be quite easy for you.
To learn how to work with anndata, you can read this `tutorial <https://falexwolf.de/blog/171223_AnnData_indexing_views_HDF5-backing/>`_ written by the author of `anndata` and `scanpy`.
You can also refer to `anndata documentation <https://anndata.readthedocs.io/en/stable/>`_.
The spatialtis analysis workflow is designed to work with `scanpy`, a python library for single cell analysis,
you can perform normal single cell analysis with `scanpy` and spatial analysis with spatialtis.
However, it's possible to skip scanpy as long as you can sort out the cell type.
If you want to use `scanpy`, please read `scanpy documentation <https://icb-scanpy.readthedocs-hosted.com/en/stable/>`_ for details,
or you can start with `scanpy tutorial <https://scanpy-tutorials.readthedocs.io/en/latest/index.html>`_.
If you are not familiar with `pandas` and `numpy`,
there are a lots of good tutorials out there, check the official tutorial for `numpy <https://numpy.org/devdocs/user/quickstart.html>`_ and `pandas <https://pandas.pydata.org/pandas-docs/stable/getting_started/tutorials.html>`_.


Reading Data
------------

spatialtis prepared a preprocessing module to help you convert your data into `anndata`.
Currently, spatialtis can read single cell information from multichannel image file with mask image prepared.
For example, technology like Imaging mass cytometry (IMC) and Multiplexed ion beam imaging (MIBI).
For other spatial single-cell technologies, it's easy to transform your expression matrix into `anndata`.

You input data must contain in *one* entry folder, each ROI must store in seperate sub-folders with its mask image,
and the structure must organized as your experiments.

Let's say your file structure look like this::

            Data
            ├── Patient1
            │   ├── Sample1
            │   │   ├── ROI1
            │   │   ├── ROI2
            │   ├──Sample2
            │   │   ├── ROI1
            │   │   ├── ROI2
            │   └──Sample3
            │       ├── ROI1
            │       └── ROI2
            ├── Patient2
            ├── Patient3
            └── ...

In one of the ROI folders, it contains data files like,
there should always be a `mask` image that tells spatialtis how each cell looks like,
spatialtis **CAN'T** do segmentation. The `mask` images should have the same naming pattern,
eg. All contain a unique keyword "mask". so that spatialtis knows which one is the `mask` image.::

    ROI1
    ├── Patient1_Body_ROI1_Dy161_CD20.tiff
    ├── Patient1_Body_ROI1_Dy162_CD8.tiff
    ├── Patient1_Body_ROI1_Dy164_CD99.tiff
    ├── Patient1_Body_ROI1_Er166_NFkB.tiff
    ├── ...
    └── Patient1_Body_ROI1_mask.tiff


Read all your data::

    entry = '/Data'
    condition_name = ['Patient', 'Sample', 'ROI']

    data = read_ROIs(entry, condition_name)

There are two type of data files, one is to store each channel in seperate images, another is to store all channels in on stacked images, if you use stacked images as input::

    data = read_ROIs(entry, condition_name, stacked=True)

To tell spatialtis the channels and markers, you can config it in two ways, you can use file::

    """metadata.csv
    channels,markers
    Dy161,CD20
    Dy162,CD8
    Dy164,CD99
    Er166,NFkB
    ...
    """
    metadata = '/metadata.csv'

    data.config_file(metadata, channel_col='channels', marker_col='markers')

Or you can directly pass in python list, but make sure the channels and markers are corresponded, and if you use stacked images, the order of channels' name will correspond to stacked order::

    channels = ['Dy161', 'Dy162', 'Dy164', 'Er166', ...]
    markers = ['CD20', 'CD8', 'CD99', 'NFkB', ...]

    data.config(channels=channels, markers=markers)

Finally, we can start processing your images into anndata::

    data = data.to_anndata()

If you have a large dataset, you can set `mp=True` to enable parallel processing, (Linux and MacOS only)::

    data = data.to_anndata(mp=True)

The default methods to determine the cell shape is "convex hull", another option is "concave hull"
(`Determine cell shape <about/implementation.html#determine-cell-shape>`_), if you want to get more accurate shape information,
you can use `polygonize="concave"`, but you need to determine the alpha value::

    data = data.to_anndata(polygonize="concave", alpha=2.0)

.. note::
    The default method

Depends on the size of your data, usually it takes 10s - 15s each ROI. Make sure to save your data on the disk::

    data.write(filename="sample.h5ad")

Let's see what's in the data::

    print(data)
    """
    AnnData object with n_obs × n_vars = 152037 × 36
        obs: 'Patient', 'Part', 'ROI', 'area', 'cell_shape', 'centroid', 'eccentricity'
        var: 'Channels', 'Markers'
    """

This means there are 152037 cells with 36 markers. In the `obs` field, 'Patients, 'Part', 'ROI' is the name for different
experiment condition, 'area', 'cell_shape', 'centroid', 'eccentricity' is calculated by spatialtis.


Basic Usage
--------------------------

Before any analysis using spatialtis, it's necessary to set up some global config, so that you don't need to specific them every time you call a function.::

    from spatialtis import CONFIG

    CONFIG.EXP_OBS = ['Patient', 'Sample', 'ROI']
    CONFIG.TYPE_COL = 'cell_type'


There are two analysis modules in spatialtis, `statistic` and `spatial`, and a visualization module `plotting`.::

    import spatialtis.sta as st
    import spatialtis.spatial as ss
    import spatialtis.plotting as sp

Now let's load the data::

    from anndata import read_h5ad

    data = read_h5ad('/sample.h5ad')

Usually an analysis function will have a corresponded visualization function, they share the same name but exists in different modules. Please don't import those function individually, it will cause conflicts.::

    # analysis function for cell components
    st.cell_components(data)
    """
    Finished!
        Add to AnnData object
        uns: 'cell_components'
    """

    # plotting function for cell components
    sp.cell_components(data)

Now that you've learn some basic of spatialtis, it can start playing around. If you want to know more about analysis and
visualization with spatialtis, go on with our `tutorial <tutorial>`_.



