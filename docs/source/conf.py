# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or mod to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import os
import sys
from pathlib import Path
import datetime
from importlib.metadata import version

import matplotlib

sys.path.insert(0, os.path.abspath("../.."))

matplotlib.use('agg')

# -- Project information -----------------------------------------------------

project = "SpatialTis"
copyright = f"{datetime.datetime.now().year}, Mr-Milk"
author = "Mr-Milk"

_current_version = '.'.join(version('spatialtis').split('.')[:3])
release = _current_version
version = _current_version

master_doc = "index"

# -- General configuration ---------------------------------------------------

# Add any Sphinx extension mod names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    "sphinx.ext.autodoc",
    # "sphinx_autodoc_typehints",
    # "sphinx.ext.napoleon",
    "sphinx.ext.autosectionlabel",
    "sphinx.ext.mathjax",
    "sphinx.ext.viewcode",
    "sphinx.ext.autosummary",
    "sphinx.ext.intersphinx",
    "numpydoc",
    "nbsphinx",
    "nbsphinx_link"
]

autosummary_generate = True

# nbsphinx_allow_errors = True
nbsphinx_input_prompt = '[%s]:'
nbsphinx_output_prompt = '[%s]:'
autodoc_typehints = "none"
autoclass_content = "class"
# autodoc_member_order = 'groupwise'
# typehints_fully_qualified = False
# Add any paths that contain templates here, relative to this directory.
templates_path = ["_templates"]
html_static_path = ['_static']
html_css_files = ['custom.css']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = []

# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
# html_theme = "sphinx_book_theme"
html_theme = "pydata_sphinx_theme"

# Add any paths that contain custom _static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin _static files,
# so a file named "default.css" will overwrite the builtin "default.css".
# html_static_path = ["_static"]
html_favicon = 'img/Logo.svg'

html_logo = 'img/Logo.svg'
# html_theme_options = {
#     "repository_url": "https://github.com/Mr-Milk/SpatialTis",
#     "use_repository_button": True,
#     "use_issues_button": True,
#     "use_edit_page_button": True,
# }

# html_theme = "#ef8992"
# html_theme_options = {
#    "rightsidebar": "true",
#    "relbarbgcolor": "#ef8992"
# }

intersphinx_mapping = {
    "seaborn": ("https://seaborn.pydata.org/", None),
    "milkviz": ("https://milkviz.readthedocs.io/en/latest/", None),
}