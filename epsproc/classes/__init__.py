
# Add some basic setup routines here...
from epsproc.util.env import isnotebook
# import xarray as xr

# Set HTML output style for Xarray in notebooks (optional), may also depend on version of Jupyter notebook or lab, or Xr
# See http://xarray.pydata.org/en/stable/generated/xarray.set_options.html
if isnotebook():
    try:
        xr.set_options(display_style = 'html')
    except NameError:
        pass
