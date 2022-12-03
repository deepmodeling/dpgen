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
import subprocess
from datetime import date
# import sys
import recommonmark
from recommonmark.transform import AutoStructify


# -- Project information -----------------------------------------------------

project = 'DP-GEN'
copyright = '2020-%d, DeepModeling' % date.today().year
author = 'DeepModeling'


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
# extensions = [
#     'recommonmark',
#     "sphinx_rtd_theme",
#     'myst_parser',
#     'sphinx_markdown_tables',
#     'sphinx.ext.autosummary'
# ]

extensions = [
    'deepmodeling_sphinx',
    'dargs.sphinx',
    "sphinx_rtd_theme",
    'myst_parser',
    'sphinx.ext.autosummary',
    'sphinx.ext.viewcode',
    'sphinxarg.ext',
    'numpydoc',
]


# Tell sphinx what the primary language being documented is.
primary_domain = 'py'

# Tell sphinx what the pygments highlight language should be.
highlight_language = 'py'

# 
myst_heading_anchors = 4

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'sphinx_rtd_theme'

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
# html_static_path = ['_static']
# html_css_files = ['css/custom.css']

autodoc_default_flags = ['members']
autosummary_generate = True
master_doc = 'index'

intersphinx_mapping = {
    "python": ("https://docs.python.org/", None),
    "dargs": ("https://docs.deepmodeling.com/projects/dargs/en/latest/", None),
    "dpdata": ("https://docs.deepmodeling.com/projects/dpdata/en/latest/", None),
    "dpdispatcher": ("https://docs.deepmodeling.com/projects/dpdispatcher/en/latest/", None),
    "ase": ("https://wiki.fysik.dtu.dk/ase/", None),
    "numpy": ("https://docs.scipy.org/doc/numpy/", None),
    "pamatgen": ("https://pymatgen.org/", None),
    "monty": ("https://guide.materialsvirtuallab.org/monty/", None),
    "paramiko": ("https://docs.paramiko.org/en/stable/", None),
    "custodian": ("https://cloudcustodian.io/docs/", None),
    "GromacsWrapper": ("https://gromacswrapper.readthedocs.io/en/latest/", None),
}


def run_apidoc(_):
    from sphinx.ext.apidoc import main
    sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
    cur_dir = os.path.abspath(os.path.dirname(__file__))
    module = os.path.join(cur_dir, "..", "dpgen")
    main(['-M', '--tocfile', 'api', '-H', 'DP-GEN API', '-o', os.path.join(cur_dir, "api"), module, '--force'])


def setup(app):
    app.connect('builder-inited', run_apidoc)
