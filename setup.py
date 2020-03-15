from setuptools import find_packages, setup

setup(name="spatialtis",
      packages=find_packages(),
      author="Mr-Milk",
      author_email="yb97643@um.edu.mo",
      python_requires='>=3.6',
      install_requires=['ray', 'anndata', 'numpy', 'pandas',
                        'scikit-image', 'pointpats', 'scipy', 'shapely',
                        'python-igraph', 'leidenalg', 'alphashape', 'bokeh',
                        'seaborn', 'colour', 'IPython', 'matplotlib', 'tqdm'])
