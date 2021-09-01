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
extensions = ['sphinx.ext.autodoc', 'sphinx.ext.napoleon',
                'sphinxcontrib.apidoc', 'recommonmark',
                'sphinx.ext.viewcode', 'nbsphinx']
                # 'IPython.sphinxext.ipython_console_highlighting']  # Actually this throws an error on RTD - try adding ipyhton to requirements.txt instead...

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
html_css_files = [
    'max_width_fix.css',
]

# This might also work:
# html_theme_options = {'body_max_width': '70%'}
