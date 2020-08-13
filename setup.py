from pathlib import Path

from setuptools import find_packages, setup

README = Path("README.md").read_text()

setup(name="spatialtis",
      packages=find_packages(include=["spatialtis"]),
      description="spatial analysis toolkit for single-cell multiplexed tissue data",
      long_description=README,
      long_description_content_type="text/markdown",
      version="1.0.0.dev0",
      author="Mr-Milk",
      url="https://github.com/Mr-Milk/SpatialTis",
      author_email="yb97643@um.edu.mo",
      license="Apache License 2.0",
      classifiers=[
          "License :: OSI Approved :: Apache Software License",
          "Programming Language :: Python :: 3",
          "Intended Audience :: Science/Research",
          "Topic :: Scientific/Engineering :: Bio-Informatics",
      ],
      python_requires='>=3.6',
      install_requires=['ray', 'anndata', 'numpy', 'pandas',
                        'scipy', 'shapely', 'bokeh',
                        'seaborn', 'colour', 'IPython', 'matplotlib', 'tqdm', 'pyecharts',
                        'spatialentropy', 'colorama', 'snapshot_phantomjs', ],
      extra_requires={'all': ['scikit-image', 'python-igraph', 'leidenalg',
                              'alphashape', 'tifffile', 'neighborhood_analysis']}
      )
