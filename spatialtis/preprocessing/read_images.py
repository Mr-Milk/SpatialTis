import warnings
from pathlib import Path
from typing import List, Optional

import pandas as pd
from anndata import AnnData
from scipy.sparse import csr_matrix

from spatialtis.typing import File
from spatialtis_core import dumps_points_wkt, dumps_polygons_wkt, points_shapes

import numpy as np

from spatialtis.utils import options_guard


def stacked_channels(c_imgs: List[File],
                     export: File = "stacked.tiff"
                     ):
    """
    A helper functions to stack single channel images to multi-channels

    Args:
        c_imgs: A list of single channel images
        export: The path to export the stacked images

    Returns:

    """
    try:
        from skimage.io import imread, imsave
    except ImportError:
        raise ImportError("Required scikit-image, try `pip install scikit-image`.")

    stacked = []
    for c in c_imgs:
        img = imread(c)
        if img.ndim != 2:
            raise ValueError(f"{c} contains more than 1 channels data")
        stacked.append(img)

    imsave(export, np.array(stacked))


def parse_img(stacked_img, img_axes: str = "cyx"):
    try:
        from tifffile import TiffFile
        from skimage.io import imread
    except ImportError:
        raise ImportError("Required scikit-image, try `pip install scikit-image`.")

    stacked_img = str(stacked_img)
    imread_status = True

    try:
        img_data = imread(stacked_img)
        if img_axes == "cyx":
            return img_data.transpose((1, 2, 0))
        return img_data

    except Exception:
        msg = f"Faild to parse {stacked_img}, it should stored as XYC, try to read from page to page"
        warnings.warn(msg)
        imread_status = False

    img_data = []
    if not imread_status:
        with TiffFile(stacked_img) as tif:
            for p in tif.pages:
                img_data.append(p.asarray())

    # TODO: 3d support?
    img_data = np.array(img_data)
    if img_data.ndim > 4:
        raise ValueError("The dimensions of image are too much")
    else:
        return img_data.transpose((1, 2, 0))


def check_exists(files):
    non_exist = []
    for f in files:
        if not Path(f).exists():
            non_exist.append(f)
    if len(non_exist) > 0:
        raise FileNotFoundError(f"{', '.join(non_exist)} not found.")


def read_ROIs(
        images: List[File],
        masks: List[File],
        markers: Optional[pd.DataFrame] = None,
        annotations: Optional[pd.DataFrame] = None,
        image_axes: str = "cyx",
        intensity_measure: str = "mean",
        shape_approx: str = "convex",
        concavity: float = 1.5,
        geopandas_compatible: bool = True,
        is3d: bool = False,
        sparse: bool = False,
) -> AnnData:

    try:
        from skimage.io import imread
        from skimage.measure import regionprops_table
    except ImportError:
        raise ImportError("Required scikit-image, try `pip install scikit-image`.")

    if isinstance(images, (Path, str)):
        images = [images]
    if isinstance(masks, (Path, str)):
        masks = [masks]

    # check file exist before continue to prevent useless work
    check_exists(images)
    check_exists(masks)

    image_axes = options_guard(image_axes, ["cyx", "xyc"])
    shape_approx = options_guard(shape_approx, ["convex", "concave"])

    properties = ['label', 'centroid', 'image_intensity',
                  'area', 'eccentricity', 'coords', 'extent',
                  'orientation', 'axis_major_length', 'axis_minor_length'
                  ]

    if annotations is None:
        annotations = pd.DataFrame({"ROI": [f"ROI{i + 1}" for i in range(len(images))]})
    roi_header = annotations.columns

    exp_data = []
    obs_data = []
    anno_collect = []

    if is3d:
        properties.remove('coords')
    for img, msk, (_, names) in zip(images, masks, annotations.iterrows()):
        img = Path(img)
        intensities = parse_img(img, img_axes=image_axes)
        msk_data = imread(msk)
        measurements = regionprops_table(msk_data, intensities, properties=properties)

        exp_all = []
        for v in measurements['image_intensity']:
            cell_exp = getattr(v, intensity_measure).__call__(axis=(0, 1))
            exp_all.append(cell_exp)
        if len(markers) != len(exp_all[0]):
            raise ValueError(f"{len(markers)} markers doesn't match {len(exp_all[0])} channels in {img}")
        meta = pd.DataFrame({
            "label": measurements['label'],
            "centroid_x": measurements['centroid-0'],
            "centroid_y": measurements['centroid-1'],
        })
        if is3d:
            meta['centroid_z'] = measurements['centroid-2']
        else:
            points = [d.values for _, d in meta[['centroid_x', 'centroid_y']].iterrows()]
            if geopandas_compatible:
                meta['wkt_centroid'] = dumps_points_wkt(points)
                shapes = points_shapes(measurements['coords'], method=shape_approx, concavity=concavity)
                meta['wkt_shape'] = dumps_polygons_wkt(shapes)

        meta['area'] = measurements['area']
        meta['major_axis'] = measurements['axis_major_length']
        meta['minor_axis'] = measurements['axis_minor_length']
        meta['eccentricity'] = measurements['eccentricity']
        meta['extent'] = measurements['extent']
        meta['orientation'] = measurements['orientation']

        # write path
        meta['source_image'] = str(img)
        meta['image_file'] = str(img.name)
        obs_data.append(meta)
        exp_data.append(exp_all)
        anno_collect += [names.values for _ in range(len(meta))]

    obs = pd.concat(obs_data).reset_index(drop=True)
    obs[roi_header] = anno_collect
    obs.index = obs.index.astype(str)
    X = np.vstack(exp_data)
    if sparse:
        X = csr_matrix(X)
    data = AnnData(obs=obs, var=markers, X=X)
    # add spatial coord to obsm
    centroid_cols = ['centroid_x', 'centroid_y', 'centroid_z'] if is3d else ['centroid_x', 'centroid_y']
    data.obsm['spatial'] = data.obs[centroid_cols].to_numpy()

    return data
