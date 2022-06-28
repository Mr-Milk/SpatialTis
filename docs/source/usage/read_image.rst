Multiplexed images to AnnData
==============================

SpatialTis could help you transform your segmented multiplexed images data
into `AnnData`. You need to prepare:

- Stacked `.tiff`-alike images
- Cell mask images

.. warning::
    SpatialTis **CAN'T** do segmentation, there are lots of resources out there
    that can help you with this part. You can try `MCMICRO <https://mcmicro.org/>`_,
    `ilastik <https://www.ilastik.org/>`_, `DeepCell <https://www.deepcell.org/>`_
    and many more.


>>> import spatialtis as st

Specify your images and masks paths. You may also specify the
marker name for each channel and add information of your images.

>>> images = ['image1.tif', 'image2.tif']
>>> masks = ['image1_mask.tif', 'image2_mask.tif']
>>> markers = pd.DataFrame(index=['CD20', 'CD8', 'CD99', 'NFkB'])
>>> annotations = pd.DataFrame({
...     "ROI": ["ROI1", "ROI2"],
...     "Tissue": ["Front", "Tail"],
... })

With all information ready, we can read it into :code:`AnnData`.
In the meantime, SpatialTis will extract the geometry information and
expression matrix from the images.

>>> data = st.read_images(images, masks, markers=markers, annotations=annotations)
>>> data
    AnnData object with n_obs × n_vars = 2037 × 4
        obs: 'Tissue', 'ROI', 'area', 'cell_shape', 'centroid', 'eccentricity', ...

