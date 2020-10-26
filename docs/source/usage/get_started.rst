Getting Started!
=================

To install SpatialTis, please see the `Installation <installation.rst>`_ section. It's recommended to start a new virtual environment before installing the SpatialTis.

Prerequisites
-------------

Must learn:

- `anndata <https://anndata.readthedocs.io/en/stable/>`_: `Tutorial <https://falexwolf.de/blog/171223_AnnData_indexing_views_HDF5-backing/>`_ written by the author of `anndata` and `scanpy`

Optional:

- `numpy <https://numpy.org/devdocs/user/quickstart.html>`_: For array operation
- `pandas <https://pandas.pydata.org/pandas-docs/stable/getting_started/tutorials.html>`_: For table operation
- `scanpy <https://scanpy-tutorials.readthedocs.io/en/latest/index.html>`_: Used for cell type calling and scRNA-seq related analysis

Basic Usage
--------------------------

Before any analysis using SpatialTis, it's necessary to set up some :doc:`global config <../api_index/config>`, so that you don't need to specific them every time you call a function::

    from spatialtis import CONFIG

    CONFIG.EXP_OBS = ['Patient', 'Sample', 'ROI']
    CONFIG.CELL_TYPE_KEY= 'cell_type'
    CONFIG.WORKING_ENV = 'jupyter_notebook'
    CONFIG.MULTI_PROCESSING = True

There are an analysis module and a visualization module `plotting` in SpatialTis::

    import spatialtis as st
    import spatialtis.plotting as sp

First load the data::

    from anndata import read_h5ad
    data = read_h5ad('/sample.h5ad')

Normally an analysis function comes with a visualization function,
they share the same name but registered in different modules.
**DON'T** import these functions individually, it will cause conflicts

Call an analysis function for cell_components::

    st.cell_components(data)

Call visualization function for cell_components::

    sp.cell_components(data)

Now that you've learn some basic of SpatialTis, feel free to play around. If you want to know more about analysis and
visualization with SpatialTis, check the :doc:`examples <example>`.

Reading Data
------------

SpatialTis prepared a preprocessing module to help you convert your data into `AnnData`.
Currently, SpatialTis can read single cell information from multichannel image file with mask image prepared.
For example, technology like Imaging mass cytometry (IMC) and Multiplexed ion beam imaging (MIBI).
For other spatial single-cell technologies, it's easy to transform your expression matrix into `AnnData`.


- Input data must contain in **one** entry folder
- Each ROI must store in separated sub-folders with its mask image
- The structure must organized as your experiments.

Let's say your file structure look like this::

            Data
            ├── Patient1
            │   ├── Sample1
            │   │   ├── ROI1
            │   │   ├── ROI2
            │   ├──Sample2
            │   │   ├── ROI1
            │   │   ├── ROI2
            │   └──Sample3
            │       ├── ROI1
            │       └── ROI2
            ├── Patient2
            ├── Patient3
            └── ...

In each one of the ROI folders contains **two** data files:

    - **One** mask image file
    - **One** stacked channels image file

There should always be a `mask` image that tells SpatialTis how each cell looks like,
SpatialTis **CAN'T** do segmentation.

If you image channels are stored in separated images, like below::

    ROI1
    ├── Patient1_Body_ROI1_Dy161_CD20.tiff
    ├── Patient1_Body_ROI1_Dy162_CD8.tiff
    ├── Patient1_Body_ROI1_Dy164_CD99.tiff
    ├── Patient1_Body_ROI1_Er166_NFkB.tiff
    ├── ...
    └── Patient1_Body_ROI1_mask.tiff

Please stacked them into single image file, your folder should eventually look
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
"channels" columns. But remember the order of your markers must align with the layers' order in your image file.
Either following the order of channels or pages, depends on the structure of your `.tiff`::

    channels = ['Dy161', 'Dy162', 'Dy164', 'Er166', ...]
    markers = ['CD20', 'CD8', 'CD99', 'NFkB', ...]
    var = pd.DataFrame({"channels":channels, "markers":markers})

Now, let's read it out::

    data = read_ROIs(entry, obs_name, var,
                     mask_pattern="mask", img_pattern="stacked")

You must noticed that there are another two arguments, *mask_pattern* allows you to tell SpatialTis which file is the mask
image, in our example, file name contains "mask" will be used as mask image. The same for *img_pattern*, so if you have
other files in the same directory, this will help SpatialTis identify which is mask image and which is data image.

Finally, we can start processing your images into anndata, the speed is related to the size of your image, in my own test
an 1000*1000 ROI from IMC data takes ~15s::

    data = data.to_anndata()

If you have a large dataset, you can set `mp=True` to enable parallel processing::

    data = data.to_anndata(mp=True)

The default methods to determine the cell shape is "convex hull", another option is "concave hull"
(:doc:`Determine cell shape <../about/implementation>`). Although "concave hull" is going to give you
more accurate shape, i strongly recommend using "convex hull".

After the processing, make sure to **save your data** on the disk::

    data.write(filename="sample.h5ad")

Let's see what's in the data::

    print(data)
    """
    AnnData object with n_obs × n_vars = 152037 × 36
        obs: 'Patient', 'Sample', 'ROI', 'area', 'cell_shape', 'centroid', 'eccentricity'
        var: 'channels', 'markers'
    """

This means there are 152037 cells with 36 markers. In the `obs` field, 'Patients, 'Sample', 'ROI' are the names for different
experiment condition, 'area', 'cell_shape', 'centroid', 'eccentricity' is calculated by SpatialTis.




