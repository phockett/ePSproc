# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# http://www.sphinx-doc.org/en/master/config

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import os
import sys
sys.path.insert(0, os.path.abspath('.'))
#sys.path.insert(0, os.path.abspath('..\\'))    # Add module path (relative to docs path)
#sys.path.insert(0, os.path.abspath('..\\..\\'))    # Add module path (relative to docs path)
sys.path.insert(0, os.path.abspath('../../'))       # Add module path (relative to docs path) FOR READTHEDOCs (above originally worked, but then broke!)
#sys.path.insert(0, os.path.abspath('..\\..\\..\\'))    # Add module path (relative to docs path)
# print(sys.path)

# Set RTD flag, use this to skip imports for RTD build.
# See https://docs.readthedocs.io/en/stable/faq.html#how-do-i-change-behavior-when-building-with-read-the-docs
on_rtd = os.environ.get('READTHEDOCS') == 'True'
# if on_rtd:
#     html_theme = 'default'
# else:
#     html_theme = 'nature'


# -- Project information -----------------------------------------------------

project = 'ePSproc'
copyright = '2019, Paul Hockett'
author = 'Paul Hockett'

# Version from package https://stackoverflow.com/questions/26141851/let-sphinx-use-version-from-setup-py
if on_rtd:
    version = "RTD"  # Dummy variable, RTD will use version from Github branch/tag.
                    # See https://readthedocs.org/projects/epsproc/versions/
else:
    from epsproc import __version__
    version = __version__

# The full version, including alpha/beta/rc tags
release = version


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
# NOTE 'IPython.sphinxext.ipython_console_highlighting' for RTD ipython highlighting.
# See https://github.com/spatialaudio/nbsphinx/issues/24
#
# Markdown support testing 09/11/22
#   OK to add ONLY ONE ext, otherwise fails on RTD, "Extension error: source_suffix '.md' is already registered".
#   Working OK with 'recommonmark' for full MD docs only, or 'sphinx_mdinclude' for full docs and rst with .. mdinclude::, see https://sphinx-mdinclude.omnilib.dev/en/latest/example.html
#   Didn't test myst_parser as yet.
#
extensions = ['sphinx.ext.autodoc', 'sphinx.ext.napoleon',
                'sphinxcontrib.apidoc', # 'recommonmark',
                'sphinx.ext.viewcode', 'nbsphinx',
                # 'myst_parser',
                'sphinx_mdinclude']    # 09/11/22 - testing MD include support, e.g. .. mdinclude:: ../../../docker/readme.md
                # 'IPython.sphinxext.ipython_console_highlighting']  # Actually this throws an error on RTD - try adding ipyhton to requirements.txt instead...

# 09/11/22 - testing MD setup, see https://github.com/readthedocs/blog/blob/main/adding-markdown-support.rst
# UPDATE: should be handled automatically by extensions now?
# source_parsers = {'.md': CommonMarkParser}
# source_suffix = ['.rst', '.md']

# api doc settings
apidoc_module_dir = '../../epsproc'
apidoc_output_dir = 'modules'
apidoc_excluded_paths = ['tests']
apidoc_separate_modules = True

# Sphinx-autodoc mock imports for minimal build-chain.
# https://www.sphinx-doc.org/en/master/usage/extensions/autodoc.html#confval-autodoc_mock_imports
if on_rtd:
    autodoc_mock_imports = ["numpy_quaternion", "quaternion", "spherical_functions","cclib",
                            "numpy","scipy","xarray","pandas","numba",
                            "matplotlib","mpl_toolkits","seaborn","plotly",
                            "pyvista","holoviews",
                            "natsort"]

    # Try fix for Plotly (currently fails on dev branch, but OK on master)
    # See https://github.com/phockett/ePSproc/issues/27
    # Fix per https://github.com/readthedocs/sphinx_rtd_theme/issues/788#issuecomment-585785027
    # Adds explict link per page
    # Appears OK in HTML, but Plotly NOT working, now added https://github.com/plotly/plotly.js/ script too...
    nbsphinx_prolog = r"""
.. raw:: html

    <script src='http://cdnjs.cloudflare.com/ajax/libs/require.js/2.1.10/require.min.js'></script>
    <script>require=requirejs;</script>
    <script src="https://cdn.plot.ly/plotly-2.24.3.min.js" charset="utf-8"></script>


"""

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ['**.ipynb_checkpoints','**-verified*','**Dev*']

# For Read the Docs, see https://stackoverflow.com/questions/56336234/build-fail-sphinx-error-contents-rst-not-found
master_doc = 'index'

# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
# html_theme = 'alabaster'
html_theme = 'sphinx_rtd_theme'

html_logo = 'figs/ePSproc_logo.png'

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']

# Add custom CSS to fix max width issue
# Fix RTD standard template (sphinx_rtd_theme) max width.
# See
# https://github.com/readthedocs/sphinx_rtd_theme/issues/295
# https://stackoverflow.com/questions/23211695/modifying-content-width-of-the-sphinx-theme-read-the-docs

# Option (1): include .css file with fix. Applied & tested on RTD OK, 01/09/21
html_css_files = [
    'max_width_fix.css',
]

# Option (2): set theme option
# This might also work (no additional CSS required):
# DOESN'T WORK ON RTD with sphinx_rtd_theme
# html_theme_options = {'body_max_width': '70%'}

# Option (3): include .css file with patch (similar to (1), but imports theme into the patch CSS file)
# This might work too - patches existing theme:
# html_style = 'max_width_patch.css'
