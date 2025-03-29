# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# from f4epurity import __version__
import sys
import os

sys.path.insert(0, os.path.abspath("../.."))

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'F4Epurity'
copyright = '2025, UKAEA - F4E'
author = 'UKAEA - F4E'
release = "v1.0.0"

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = ["myst_parser", "sphinx.ext.autodoc", "sphinx.ext.napoleon"]

source_suffix = {
    ".rst": "restructuredtext",
    ".md": "markdown",
}
# Add any paths that contain templates here, relative to this directory.
templates_path = ["_templates"]

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = []

# Option for autodoc
autodoc_default_options = {"show-inheritance": True, "exclude-members": "__weakref__"}
autoclass_content = "init"

pygments_style = "sphinx"
# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = "pydata_sphinx_theme"
#html_logo = "png_path"
# html_theme_options = {
#    "github_url": "repo_url",
#}