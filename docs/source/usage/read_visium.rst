Read Visium Result Folders to AnnData
======================================

After processing your visium experiments, you should have
a folder for each ROI. You can read all of them into
one `AnnData`:

>>> import spatialtis as st
>>> folders = [
...   "ROI_1_filtered_gene_bc_matrix_spatial",
...   "ROI_2_filtered_gene_bc_matrix_spatial",
...   ]
>>> annotations = pd.DataFrame({
...     "ROI": ["ROI1", "ROI2"],
...     "Tissue": ["Front", "Tail"],
... })

By default it will read from the filtered matrix

>>> data = st.read_visium(folders, read_filtered=True, annotations=annotations)
>>> data
    AnnData object with n_obs × n_vars = 8072 × 10467
        obs: 'Tissue', 'ROI'
        obsm: 'spatial'

