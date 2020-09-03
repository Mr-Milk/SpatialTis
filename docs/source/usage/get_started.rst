Getting Started!
=================

To install spatialtis, please see the `Installation <installation.rst>`_ section. It's recommended to start a new virtual environment before installing the spatialtis.

Prerequisites
-------------

To work with spatialtis, you need to know about `anndata`, if you are familiar with `pandas` and `numpy` this should be quite easy for you.
To learn how to work with anndata, you can read this `tutorial <https://falexwolf.de/blog/171223_AnnData_indexing_views_HDF5-backing/>`_ written by the author of `anndata` and `scanpy`.
You can also refer to `anndata documentation <https://anndata.readthedocs.io/en/stable/>`_.

If you are not familiar with `pandas` and `numpy`,
there are a lots of good tutorials out there, check the official tutorial
for `numpy <https://numpy.org/devdocs/user/quickstart.html>`_
and `pandas <https://pandas.pydata.org/pandas-docs/stable/getting_started/tutorials.html>`_.

The spatialtis analysis workflow is designed to work with `scanpy`, a python library for single cell analysis,
you can perform normal single cell analysis with `scanpy` and spatial analysis with spatialtis.
However, it's possible to skip scanpy as long as you can sort out the cell type.
If you want to use `scanpy`, please read `scanpy documentation <https://icb-scanpy.readthedocs-hosted.com/en/stable/>`_ for details,
or you can start with `scanpy tutorial <https://scanpy-tutorials.readthedocs.io/en/latest/index.html>`_.

Basic Usage
--------------------------

Before any analysis using spatialtis, it's necessary to set up some :doc:`global config <../api_index/config>`, so that you don't need to specific them every time you call a function.::

    from spatialtis import CONFIG

    CONFIG.EXP_OBS = ['Patient', 'Sample', 'ROI']
    CONFIG.TYPE_COL = 'cell_type'


There are two analysis modules in spatialtis, `statistic` and `spatial`, and a visualization module `plotting`.::

    import spatialtis as st
    import spatialtis.plotting as sp

Now let's load the data::

    from anndata import read_h5ad
    data = read_h5ad('/sample.h5ad')

Usually an analysis function will have a corresponded visualization function,
they share the same name but exists in different modules.
Please don't import those function individually, it will cause conflicts.::

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
visualization with spatialtis, go on with our tutorial.

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

In each one of the ROI folders, show contains *two* data files:

    - *One* mask image file
    - *One* stacked channels image file

there should always be a `mask` image that tells spatialtis how each cell looks like,
spatialtis **CAN'T** do segmentation.

If you image channels are stored in seperated images, like below::

    ROI1
    ├── Patient1_Body_ROI1_Dy161_CD20.tiff
    ├── Patient1_Body_ROI1_Dy162_CD8.tiff
    ├── Patient1_Body_ROI1_Dy164_CD99.tiff
    ├── Patient1_Body_ROI1_Er166_NFkB.tiff
    ├── ...
    └── Patient1_Body_ROI1_mask.tiff

please stacked them into single image file, your folder should eventually look
like this::

    ROI1
        ├── stacked.tiff
        ├── Patient1_Body_ROI1_mask.tiff
        └── ...

First we need to specific the path of the entry folder of our data::

    entry = '/Data'

And then, we need to describe how the experiment is designed.
The names should be corresponded to each level of the folder, if you look back to the file tree
that we show before, it's easy to understand::

    obs_name = ['Patient', 'Sample', 'ROI']

Another information is the markers data, it needs to be stored into a dataframe.
This allow you to add as many columns as you want, for example you can add an extra
"channels" columns. But remember the order of your marker must correspond to the layers in your image file.
Either following the order of channels or pages, depends on the structure of your `.tiff`::

    channels = ['Dy161', 'Dy162', 'Dy164', 'Er166', ...]
    markers = ['CD20', 'CD8', 'CD99', 'NFkB', ...]
    var = pd.DataFrame({"channels":channels, "markers":markers})

Now, let's read it out::

    data = read_ROIs(entry, obs_name, var,
                     mask_pattern="mask", img_pattern="stacked")

You must noticed that there are another two arguments, *mask_pattern* allows you to tell spatialtis which file is the mask
image, in our example, file name contains "mask" will be used as mask image. The same for *img_pattern*, so if you have
other files in the same directory, this will help spatialtis identify which is mask image and which is data image.

Finally, we can start processing your images into anndata, the speed is related to the size of your image, in my own test
an 1000*1000 ROI from IMC data takes ~15s::

    data = data.to_anndata()

If you have a large dataset, you can set `mp=True` to enable parallel processing::

    data = data.to_anndata(mp=True)

The default methods to determine the cell shape is "convex hull", another option is "concave hull"
(:doc:`Determine cell shape <../about/implementation>`). Although "concave hull" is going to give you
more accurate shape, i strongly recommend using "convex hull".

After the processing Make sure to save your data on the disk::

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




