Multiplexed images to AnnData
==============================

For image-based spatial single-cell technologies, SpatialTis could help you transform
your images data into `AnnData`. You need to prepare:

- Stacked `.tiff` images
- Cell mask image

.. warning::
    SpatialTis **CAN'T** do segmentation

The file structure should look like::

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

The minimum folder is the ROI folder, each contains **two** data files::

    ROI1
    ├── stacked.tiff
    └── Patient1_Body_ROI1_mask.tiff

- **One** stacked image
- **One** mask image

>>> import spatialtis as st

First, we need to specify the entry points

>>> entry = "./Data"

And then, we need to describe how the experiment is designed.
The names should be corresponded to each level of the folder,
if you look back to the file tree that we show before, it’s easy to understand

>>> obs_name = ['Patient', 'Sample', 'ROI']

Another information is the markers data,
it needs to be stored into a dataframe.
This allow you to add as many columns as you want,
for example you can add an extra “channels” columns.
But remember the order of your markers must align with the
layers’ order in your image file.
Either following the order of channels or pages, depends on the structure of your *.tiff*:

>>> channels = ['Dy161', 'Dy162', 'Dy164', 'Er166', ...]
>>> markers = ['CD20', 'CD8', 'CD99', 'NFkB', ...]
>>> var = pd.DataFrame({"channels":channels, "markers":markers})

Now, let’s read it out

>>> reader = read_ROIs(entry, obs_name, var, mask_pattern="mask", img_pattern="stacked")

`mask_pattern` allows you to tell SpatialTis which file is the mask image,
in our example, file name contains “mask” will be used as mask image.
The same for `img_pattern`, so if you have other files in the same directory,
this will help SpatialTis identify which is mask image and which is data image.

Finally, we can start processing your images into anndata,
the speed is related to the size of your image,
in my own test an 1000*1000 ROI from IMC data takes ~15s.

>>> data = reader.to_anndata()
>>> data.write("sample.h5ad")


>>> print(data)
    AnnData object with n_obs × n_vars = 152037 × 36
        obs: 'Patient', 'Sample', 'ROI', 'area', 'cell_shape', 'centroid', 'eccentricity'
        var: 'channels', 'markers'

