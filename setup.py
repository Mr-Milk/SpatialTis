from setuptools import setup, find_packages

setup(name="spatialtis", packages=find_packages(), install_requires=['ray', 'anndata', 'numpy', 'pandas',
                                                                     'scikit-image', 'pointpats', 'scipy', 'shapely',
                                                                     'igraph', 'leidenalg', 'alphashape', 'bokeh',
                                                                     'seaborn', 'colour'])
