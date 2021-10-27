# import os
# from pathlib import Path
#
# from setuptools import find_packages, setup
# import versioneer
#
#
# README = Path("README.md").read_text(encoding='utf-8')
#
# setup(name="spatialtis",
#       packages=find_packages(),
#       description="spatial analysis toolkit for single-cell multiplexed tissue data",
#       long_description=README,
#       long_description_content_type="text/markdown",
#       version=versioneer.get_version(),
#       author="Mr-Milk",
#       url="https://github.com/Mr-Milk/SpatialTis",
#       author_email="yb97643@um.edu.mo",
#       license="Apache License 2.0",
#       classifiers=[
#           "License :: OSI Approved :: Apache Software License",
#           "Programming Language :: Python :: 3",
#           "Programming Language :: Python :: 3.6",
#           "Programming Language :: Python :: 3.7",
#           "Programming Language :: Python :: 3.8",
#           "Programming Language :: Python :: 3.9",
#           "Intended Audience :: Science/Research",
#           "Topic :: Scientific/Engineering :: Bio-Informatics",
#       ],
#       python_requires='>=3.6',
#       install_requires=['anndata', 'numpy', 'pandas', 'scipy', 'scikit-learn', 'rich[jupyter]',
#                         'spatialtis_core', 'milkviz'],
#       extras_require={'all': ['scikit-image', 'ray', 'lightgbm',
#                               'leidenalg', 'python-igraph', ]}
#       )

import setuptools

setuptools.setup()
