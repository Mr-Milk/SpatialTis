# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import os
import sys

sys.path.insert(0, os.path.abspath("../.."))

# -- Project information -----------------------------------------------------

project = "SpatialTis"
copyright = "2020, Mr-Milk"
author = "Mr-Milk"

# The full version, including alpha/beta/rc tags
release = "0.3.0"

master_doc = "index"

# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    "sphinx.ext.autodoc",
    "sphinx_autodoc_typehints",
    "sphinx.ext.napoleon",
    "sphinx.ext.autosectionlabel",
    "sphinx.ext.mathjax",
    "sphinx.ext.viewcode",
    "sphinx_rtd_theme",
    "nbsphinx",
    "nbsphinx_link"
]

# nbsphinx_allow_errors = True
# nbsphinx_input_prompt = ''
# nbsphinx_output_prompt = ''
autodoc_typehints = "description"
autoclass_content = "both"
autodoc_member_order = 'groupwise'
typehints_fully_qualified = False
# Add any paths that contain templates here, relative to this directory.
templates_path = ["_templates"]
html_static_path = ['_static']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = []

# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = "sphinx_rtd_theme"

# Add any paths that contain custom _static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin _static files,
# so a file named "default.css" will overwrite the builtin "default.css".
# html_static_path = ["_static"]
html_favicon = 'img/Logo.svg'

html_logo = 'img/Logo-text.png'
html_theme_options = {
    'logo_only': True,
}

# html_theme = "#ef8992"
# html_theme_options = {
#    "rightsidebar": "true",
#    "relbarbgcolor": "#ef8992"
# }
