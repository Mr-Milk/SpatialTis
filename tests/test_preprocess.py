import pandas as pd
import spatialtis as st


def test_read_images(ome_images):
    images, masks = ome_images
    annotations = pd.DataFrame({
        "ROI": ["ROI1", "ROI2"]
    })
    st.read_images(images, masks, sparse=True, annotations=annotations)
    st.read_images(images, masks, shape_approx="concave")


def test_read_visium(visium_folder):
    st.read_visium(visium_folder)
