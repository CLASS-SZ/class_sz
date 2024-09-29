# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'class_sz'
copyright = '2024, Bolliet, Kusiak, McCarthy et al'
author = 'Bolliet, Kusiak, McCarthy et al'


import os
import tomli


# Extract version from pyproject.toml
def get_version_from_pyproject():
    with open(os.path.join(os.path.dirname(__file__), '..', 'pyproject.toml'), 'rb') as f:
        data = tomli.load(f)
        return data['project']['version']

try:
    release = get_version_from_pyproject()
except (FileNotFoundError, KeyError):
    release = '0.0.0'  # Fallback version

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration


# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
	'nbsphinx',
	'sphinx.ext.mathjax',
	'sphinx_rtd_theme',
    'sphinx_gallery.load_style',  # load CSS for gallery (needs SG >= 0.6)
    'myst_parser'
]

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']



# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']


source_suffix = ['.rst', '.md']
