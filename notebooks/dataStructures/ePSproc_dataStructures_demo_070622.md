---
jupytext:
  formats: ipynb,md:myst
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
    jupytext_version: 1.13.7
kernelspec:
  display_name: Python 3 (ipykernel)
  language: python
  name: python3
---

# Data stuctures - basic overview and demo
07/06/22

This notebook extends the [basic overview](https://epsproc.readthedocs.io/en/dev/demos/ePSproc_demo_Aug2019.html#Structure), including some updated functionality and file IO.

Note that the focus here is on low-level functions and base data structure handling, see the [class demo](https://epsproc.readthedocs.io/en/dev/demos/ePSproc_class_demo_161020.html) for more general usage, and class data structures.

Firstly, load some demo data to play with...

+++

## Load data

```{code-cell} ipython3
from pathlib import Path

import epsproc as ep

# Set data path
# Note this is set here from ep.__path__, but may not be correct in all cases - depends on where the Github repo is.
epDemoDataPath = Path(ep.__path__[0]).parent/'data'
```

```{code-cell} ipython3
# Load data from modPath\data
dataPath = Path(epDemoDataPath, 'photoionization')
dataFile = Path(dataPath, 'n2_3sg_0.1-50.1eV_A2.inp.out')  # Set for sample N2 data for testing

# Scan data file
dataSet = ep.readMatEle(fileIn = dataFile.as_posix())
data = dataSet[0]
# dataXS = ep.readMatEle(fileIn = dataFile.as_posix(), recordType = 'CrossSection')  # XS info currently not set in NO2 sample file.
```

## Xarray for ND data handling

All the core data is handled as [Xarrays](http://xarray.pydata.org/en/stable/index.html). This provides a wrapped for numpy ND arrays, including labelled coordinates, and various computational functions. For more details, see the [Xarray Data Structures documentation](https://docs.xarray.dev/en/latest/user-guide/data-structures.html).

```{code-cell} ipython3
# Calling the array will provide some readable summary output
data
```

### Xarray functionality

Various low-level functions are available...

+++

Subselect data, see https://docs.xarray.dev/en/latest/user-guide/indexing.html#indexing-and-selecting-data

```{code-cell} ipython3
inds = {'Type':'L','Cont':'PU','mu':1}  # Set a dictionary of indexes (dimensions & coordinate labels) for selection
data.sel(inds).squeeze(drop=True).dropna(dim='LM',how='all')  # Select & drop redundant coords
```

Standard max, min etc. functionality...

```{code-cell} ipython3
data.max()
```

```{code-cell} ipython3
data.min()
```

```{code-cell} ipython3
data.mean()
```

... and dimension labels can be used here ...

```{code-cell} ipython3
# Get max value over Eke dim, then sum over other dims.
# data.max(dim='Eke').sum(['Type','Sym','mu'])

# Sum over some dims, then drop any singleton dims
data.sum(dim=['Type','Sym','mu']).squeeze()
```

Data is available as a numpy ND array at `.values`

```{code-cell} ipython3
type(data.values)
```

As well as the core Xarray functionality, data can be piped directly to [any numpy universal function](https://docs.xarray.dev/en/latest/user-guide/computation.html) ...

```{code-cell} ipython3
import numpy as np
data.sum(dim=['Type','Sym','mu']).squeeze().pipe(np.abs)
```

#### Coordinates

Coordinates are also Xarrays, and accessible via dot notation for individual items, or for all at `.coords`.

See the [Xarray coords documentation](https://docs.xarray.dev/en/latest/user-guide/data-structures.html#coordinates) for more details.

Note that multiindex coordinates are also supported, as are "non-dimensional" coordinates, which can be used to provide alternative labels for existing dimensions.

```{code-cell} ipython3
# Single index numerical coordinate
data.Eke
```

```{code-cell} ipython3
# Multiindex coordinate
data.Sym
```

```{code-cell} ipython3
# All coordinates.
# Here multidimensional coords are `-` and non-dimensional coordinates are unmarked.
data.coords
```

```{code-cell} ipython3
# For low-level manipulation, these can be returned Pandas index objects.
data.Eke.to_index()
```

```{code-cell} ipython3
data.Sym.to_index()
```

### Wrapped functionality

For specific functionality, ePSproc has higher-level wrapper for various data-manuipulation tasks built on Xarray, Numpy and Pandas.

```{code-cell} ipython3
# Matrix element selector wraps Xarray selection routines & thresholding
ep.matEleSelector(data, thres=1e-2, inds = inds, sq = True)
```

For quick tabulation, any ND array can be restacked and pushed to a Pandas DataFrame.

```{code-cell} ipython3
dataRed = ep.matEleSelector(data, thres=1e-1, sq = True)  # Threshold
dataPD, _ = ep.multiDimXrToPD(dataRed, colDims = 'Eke')  # Convert to PD
dataPD
```

## Data models & types

Some core data types are defined in `ep.util.listFuncs.dataTypesList()`. These provide a reference for other functionality - although, in general, the use of Xarrays should make most routines agnostic to dimension names and ordering, some routines look for specific dimensions.

```{code-cell} ipython3
# Calling directly returns a full dictionary, or a specific dataType can be requested by key
ep.util.listFuncs.dataTypesList()['matE']
```

```{code-cell} ipython3
# The 'def' field references the base dataType definition for more options, e.g. to check stacked (multiindex) vs. unstacked definitions
print(ep.util.listFuncs.dataTypesList()['matE']['def'](sType = 'stacked'))
print(ep.util.listFuncs.dataTypesList()['matE']['def'](sType = 'unstacked'))
```

```{code-cell} ipython3
# The `sDict` type returns the stacked dim mappings
ep.util.listFuncs.dataTypesList()['matE']['def'](sType = 'sDict')
```

Based on this, it is easy to swap and rearrage dimensions (again, see the [Xarray documentation for more details](https://docs.xarray.dev/en/latest/user-guide/reshaping.html)), and some wrappers are provided for higher-level use.

```{code-cell} ipython3
# Unstack all dims with Xarray functionality
dataUnstacked = data.unstack()
dataUnstacked.coords
```

The `ep.util.misc.restack()` routine wraps some core functionality to restack arrays according to defined ePSproc dataTypes, and will try and restack according to `self.attrs['dataType']` by default.

```{code-cell} ipython3
# "Safe" dimension restacker
# This will always skip missing dims rather than throwing errors.
dataRestacked, dims = ep.util.misc.restack(dataUnstacked)
dataRestacked.coords
```

```{code-cell} ipython3
# Try restacking to a different dataType
# The default behaviour here is to skip/ignore missing dimensions
dataRestacked, dims = ep.util.misc.restack(dataUnstacked, refDims = 'BLM')
dataRestacked.coords
```

```{code-cell} ipython3
# The returned dims variable is a dictionary with various dim lists, differences and intersections, which is used by the restacked
# TODO: naming and definitions here!
# See ep.util.misc.checkDims()
dims
```

```{code-cell} ipython3
# Restack with extra stacked dims - forceUnstack = False
# This will leave extra stacked dims intact
# Note this will fail if the "extra" stacked dims include dims which should be stacked elsewhere according to refDims.
daR, _ = ep.util.misc.restack(data.stack({'Test':['Eke','it']}), forceUnstack = False)  # OK!
# daR, _ = ep.util.misc.restack(dataTest.stack({'Test':['Eke','it','Targ']}), dataTypesList()[dataType]['def'](sType='sDict'), forceUnstack = False)  # FAILS - 'Targ' in both 'Test' stack, and refDims.
                                                                                                                                                    # This case is currently not defined in checkDims()
daR.coords
```

```{code-cell} ipython3
# Restack with missing dims

# Copy data & drop a dimension
dataTest = data.copy()
dataTest['Sym'] = dataTest.indexes['Sym'].droplevel('Cont')

daR, _ = ep.util.misc.restack(dataTest)  # With missing dims
daR.coords
```

```{code-cell} ipython3
# Restack with missing dims & add dims
# Note added dims have unset coord values, but will be correctly dimensioned
daR, _ = ep.util.misc.restack(dataTest, conformDims = True)  # With missing dims
daR.coords
```

+++ {"tags": []}

## Basic data IO (Xarray data file read/write)

For Xarrays, these can be written to disk and read back with `ep.IO.writeXarray()` and `ep.IO.readXarray()`. These make use of [Xarray's netCDF writer](https://docs.xarray.dev/en/latest/user-guide/io.html#netcdf), and wrap some additional complex number and dim handling for ePSproc cases. 

For other low-level data IO options, see the [Xarray documentation](https://docs.xarray.dev/en/latest/user-guide/io.html). Higher level routines are currently in development for ePSproc & PEMtk, see further notes on [ePSproc issues](https://github.com/phockett/ePSproc/issues/8) and [PEMtk issues](https://github.com/phockett/PEMtk/issues/6) for ongoing issues and changes.

UPDATE 30/06/22: still some work to do on IO, but [more file IO details and options are now implemented and discussed on the IO page](https://epsproc.readthedocs.io/en/dev/dataStructures/ePSproc_dataStructures_IO_demo_280622.html).

```{code-cell} ipython3
dataPath = Path(epDemoDataPath, 'photoionization')
dataFile = Path(dataPath, 'n2_3sg_0.1-50.1eV_A2.nc')  # Set for sample N2 data for testing

ep.IO.writeXarray(data, fileName = dataFile.as_posix(), filePath = dataPath.as_posix())   # Default case set as: engine = 'h5netcdf', forceComplex = False
```

```{code-cell} ipython3
dataIn = ep.IO.readXarray(fileName = dataFile.as_posix() + '.nc', filePath = dataPath.as_posix())  #, forceComplex=forceComplex, forceArray=False)
```

```{code-cell} ipython3
dataIn
```

```{code-cell} ipython3
# Testing for equality with the original currently returns false - dim ordering might be the issue here?
dataIn.equals(data)
```

```{code-cell} ipython3
# Subtraction indicates identical data however.
(dataIn - data).max()
```

### Low-level routines

The Xarray routines can be used directly, although may require additional args. Note that the netcdf writer always writes to DataSet format, although there are readers for both dataset and dataarrays.

```{code-cell} ipython3
# Read file with base open_dataset
import xarray as xr
dataIn = xr.open_dataset(dataFile.as_posix() + '.nc', engine = 'h5netcdf')
dataIn
```

```{code-cell} ipython3
# Read file with base open_dataarray - this may fail for datasets
dataIn = xr.open_dataarray(dataFile.as_posix() + '.nc', engine = 'h5netcdf')
dataIn
```

### Complex number and attribute handling

Xarray methods for complex data have some limitations...

1. Issues with complex data, which is not supported by netcdf - either convert to Re+Im format or use h5netcdf backend with `forceComplex=True` to workaround (but need to set this on file read too).
1. Issues with nested dict attribs.
1. Issues with tuples, esp. in Euler dim coords.

With h5netcdf and `invalid_netcdf` option (1) is OK (see [Xarray docs for details](https://docs.xarray.dev/en/latest/user-guide/io.html#invalid-netcdf-files)), although still needs unstack, and may also need to set 'to_dataset' for more control, otherwise can get arb named items in file (if dataarray name is missing).

For general attrib handling, `ep.IO.sanitizeAttrsNetCDF()` attempts to quickly clear this up at file IO if there is an exception raised, although may be lossy in some cases.

TODO: 

- further work on this, and alternative file writers, see notes at https://github.com/phockett/ePSproc/issues/8 and https://github.com/phockett/PEMtk/issues/6.
- more attrib testing.


+++

## Versions

```{code-cell} ipython3
import scooby
scooby.Report(additional=['epsproc', 'holoviews', 'hvplot', 'xarray', 'matplotlib', 'bokeh'])
```

```{code-cell} ipython3
# Check current Git commit for local ePSproc version
from pathlib import Path
!git -C {Path(ep.__file__).parent} branch
!git -C {Path(ep.__file__).parent} log --format="%H" -n 1
```

```{code-cell} ipython3
# Check current remote commits
!git ls-remote --heads git://github.com/phockett/ePSproc
```

d87199802bc7f64cde181fedb08a614be9de6a24
