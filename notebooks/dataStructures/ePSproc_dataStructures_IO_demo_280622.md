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

# Data structures - IO demo

28/06/22

This notebook outlines various IO methods for reading and writing data. For a general datastructures overview, see [the Data structures intro doc](https://epsproc.readthedocs.io/en/latest/dataStructures/ePSproc_dataStructures_demo_070622.html).

Note that this page details functional forms, for class usage see [the base class intro page](https://epsproc.readthedocs.io/en/latest/demos/ePSproc_class_demo_161020.html). However, as of June 2022, file writers are not fully implemented for the data class.

+++

## Load libraries

```{code-cell} ipython3
from pathlib import Path

import epsproc as ep

# Set data path
# Note this is set here from ep.__path__, but may not be correct in all cases - depends on where the Github repo is.
epDemoDataPath = Path(ep.__path__[0]).parent/'data'
```

## Loading ePolyScat data

To start with, load data from an ePolyScat output file. This will load matrix elements & cross-sections, for more details see [the basic demo page](https://epsproc.readthedocs.io/en/latest/demos/ePSproc_demo_Aug2019.html) and [the ePolyScat basics page](https://epsproc.readthedocs.io/en/latest/ePS_ePSproc_tutorial/ePS_tutorial_080520.html#Results).

```{code-cell} ipython3
# Load data from modPath\data
dataPath = Path(epDemoDataPath, 'photoionization')
dataFile = Path(dataPath, 'n2_3sg_0.1-50.1eV_A2.inp.out')  # Set for sample N2 data for testing

# Scan data file
dataSet = ep.readMatEle(fileIn = dataFile.as_posix())
data = dataSet[0]
# dataXS = ep.readMatEle(fileIn = dataFile.as_posix(), recordType = 'CrossSection')  # XS info currently not set in NO2 sample file.
```

```{code-cell} ipython3
# All data is pushed to Xarrays
data
```

+++ {"tags": []}

## Xarray methods

For Xarrays, these can be written to disk and read back with `ep.IO.writeXarray()` and `ep.IO.readXarray()`. 

These wrap various backends, including [Xarray's netCDF writer](https://docs.xarray.dev/en/latest/user-guide/io.html#netcdf), and implement some additional complex number and dim handling for ePSproc cases (depending on the backend). In general, the file writers aim to preserve all dimensions, including stacking, if possible.

For other low-level data IO options, see the [Xarray documentation](https://docs.xarray.dev/en/latest/user-guide/io.html). Higher level routines are currently in development for ePSproc & PEMtk, see further notes at https://github.com/phockett/ePSproc/issues/8 and https://github.com/phockett/PEMtk/issues/6.

+++ {"jp-MarkdownHeadingCollapsed": true, "tags": []}

### Complex number and attribute handling

Note that native Xarray methods have some limitations...

1. Issues with complex data, which is not supported by netcdf - either convert to Re+Im format or use h5netcdf backend with `forceComplex=True` to workaround (but need to set this on file read too).
1. Issues with nested dict attribs.
1. Issues with MultiIndex coords, esp. in non-dim coords.

With h5netcdf and `invalid_netcdf` option (1) is OK (see [Xarray docs for details](https://docs.xarray.dev/en/latest/user-guide/io.html#invalid-netcdf-files)), although still needs unstack, and may also need to set 'to_dataset' for more control, otherwise can get arb named items in file (if dataarray name is missing).

For general attrib handling, `ep.IO.sanitizeAttrsNetCDF()` attempts to quickly clear this up at file IO if there is an exception raised, although may be lossy in some cases.

For non-Xarray file types, e.g. HDF5, some additional custom handling is included, mainly via dictionary conversion routines, which are discussed below.

TODO: 

- further work on this, and alternative file writers, see notes at https://github.com/phockett/ePSproc/issues/8 and https://github.com/phockett/PEMtk/issues/6.
- more attrib testing.


+++

### NetCDF methods

For [Xarray's netCDF writer](https://docs.xarray.dev/en/latest/user-guide/io.html#netcdf), valid engines (handlers/libraries/formats) are `h5netcdf` (preferred), `netcdf4` and `scipy`. The default case uses `engine = 'h5netcdf'`, and `forceComplex = False` - the output file will include separate real and imaginary components and a flat representation of the array.

TODO: test for multiple Xarray per file.

```{code-cell} ipython3
dataPath = Path(epDemoDataPath, 'photoionization','fileIOtests')
dataFile = Path(dataPath, 'n2_3sg_0.1-50.1eV_A2')  # Set for sample N2 data for testing

ep.IO.writeXarray(data, fileName = dataFile.as_posix(), filePath = dataPath.as_posix())   # Default case set as: engine = 'h5netcdf', forceComplex = False
```

```{code-cell} ipython3
dataIn = ep.IO.readXarray(fileName = dataFile.as_posix() + '.nc', filePath = dataPath.as_posix())  #, forceComplex=forceComplex, forceArray=False)
```

```{code-cell} ipython3
dataIn  # Note dims are restacked, although ordering may change
```

```{code-cell} ipython3
# Testing for equality with the original currently returns false - dim ordering might be the issue here?
dataIn.equals(data)
```

```{code-cell} ipython3
# Subtraction indicates identical data however.
(dataIn - data).max()
```

```{code-cell} ipython3
# For the 
```

+++ {"tags": []}

### Low-level Xarray routines

The [Xarray native routines](https://docs.xarray.dev/en/latest/user-guide/io.html) can be used directly, although may require additional args and/or data reformatting.

Note that the netcdf writer always writes to DataSet format, although there are readers for both dataset and dataarrays. Other available methods include Pickle, Zarr...

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

### HDF5

HDF5 support is available if the `h5py` library is present, and is wrapped for Xarray use with the routines in `epsproc.ioBackends.hdf5IO`. These handle most python native and numpy data types. The routines also include string wrappers for unsupported data types (and/or they can be purged for saving). For more details see the dictionary conversion options section below (the Xarray objects are serialized to dictionaries prior to file IO).

The main benefit of HDF5 (vs. netCDF) is more type handling (esp. complex data), and cross-compatibility with other languages. String wrapping is currently applied to all fields except the numerical data, ensuring handling of nested dictionaries. (See source for [writeXarrayToHDF5](https://github.com/phockett/ePSproc/blob/184abf8f06e4562b3d582dac5bd6931cef7321dc/epsproc/ioBackends/hdf5IO.py#L67) and [sanitizeAttrsNetCDF](https://github.com/phockett/ePSproc/blob/184abf8f06e4562b3d582dac5bd6931cef7321dc/epsproc/util/xrIO.py#L36).)

NOTE: string unwrapping is done via [ast.literal_eval()](https://docs.python.org/3.8/library/ast.html#ast.literal_eval), which supports: strings, bytes, numbers, tuples, lists, dicts, sets, booleans, and None. This should be safe in general (cannot evaluate arbitrary code), but pass `evalStrings = True` at file read to avoid if required.

TODO: 

- more options here, nested dict unpacking to HDF5 as well as string wrap (see .
- fix XR name, appears as byte string
- Test multiple items per file.
- Test IO with, e.g., Matlab.


Other libs/tools for this to be investigated... 

- General discussion of dicts to HDF5: https://stackoverflow.com/questions/16494669/how-to-store-dictionary-in-hdf5-dataset
   - Native methods: serialise to JSON, wrap as string.
   - Either method should be easy to implement.
- Maybe Benedict for dict handling? https://github.com/fabiocaccamo/python-benedict#io-methods
- HDFdict for nested dicts to HDF5: https://pypi.org/project/hdfdict/  (small, dicts only)
- Deepdish: https://deepdish.readthedocs.io/en/latest/io.html (general)

```{code-cell} ipython3
dataPath = Path(epDemoDataPath, 'photoionization','fileIOtests')
dataFile = Path(dataPath, 'n2_3sg_0.1-50.1eV_A2.h5')  # Set for sample N2 data for testing

ep.IO.writeXarray(data, fileName = dataFile.as_posix(), filePath = dataPath.as_posix(), engine = 'hdf5')   # Default case set as: engine = 'h5netcdf', forceComplex = False
```

```{code-cell} ipython3
dataPath = Path(epDemoDataPath, 'photoionization','fileIOtests')
dataFile = Path(dataPath, 'n2_3sg_0.1-50.1eV_A2.h5')  # Set for sample N2 data for testing

# HDF5 file reader, note this returns both dict and Xarray formats
dataInDict, dataInXR = ep.IO.readXarray(fileName = dataFile.as_posix(), filePath = dataPath.as_posix(), engine = 'hdf5')
```

```{code-cell} ipython3
# Reconstructed array. Note this also now includes self.attrs['dimMaps'], which was added by the file writer to track dimensions.
dataInXR
```

```{code-cell} ipython3
# data.attrs['dimMaps'] = dataInXR.attrs['dimMaps']
print(dataInXR.equals(data))  # False...
(dataInXR - data).max()   # But data seems OK
```

```{code-cell} ipython3
# Additional data can be appended to the file
# Note that the top-level group defaults to data.dataType, so a different name may be required for this
ep.IO.writeXarray(data, fileName = dataFile.as_posix(), filePath = dataPath.as_posix(), engine = 'hdf5', dataName = 'matE2')
```

### Low-level HDF5 IO

Native HDF5 file IO can also be used (and the files should be readable by any standard HDF5 tools), and this contains a dictionary-like representation of the data structure. For more details on [h5py see the docs](https://docs.h5py.org/en/stable/quick.html).

```{code-cell} ipython3
import h5py

dataPath = Path(epDemoDataPath, 'photoionization','fileIOtests')
dataFile = Path(dataPath, 'n2_3sg_0.1-50.1eV_A2.h5') 

# Open file with h5py
fIn = Path(dataPath,dataFile)
hf = h5py.File(fIn.as_posix(), 'r')
```

```{code-cell} ipython3
# Top-level groups
hf.keys()
```

```{code-cell} ipython3
# Group items
hf['matE'].keys()
```

```{code-cell} ipython3
# Get data from a subgroup
hf['matE']['name'][()]
```

```{code-cell} ipython3
# Get data from a subgroup
hf['matE']['dims'][()]
```

```{code-cell} ipython3
# Close file stream
hf.close()
```

## Numpy and Pandas

Conversion to other classes can be used to access their native IO routines. Note, however, that this may be lossy, depending on the type of object converted.

+++

### Pandas

Convert to a Pandas DataFrame with `ep.multiDimXrToPD()`, then use [any standard Pandas IO functionality](https://pandas.pydata.org/docs/user_guide/io.html).

In general this shoud work for all cases (including complex data), but may lose metadata/attributes.

For general data sharing, `to_hdf` and `to_csv` methods are quite useful. Note Pandas `to_hdf` [uses PyTables on the backend](https://pandas.pydata.org/docs/user_guide/io.html#hdf5-pytables).

TODO: wrapper for attributes & restacking? See PEMtk methods.

```{code-cell} ipython3
# Convert to Pandas dataframe & save to HDF5
dataFile = Path(dataPath, 'n2_3sg_0.1-50.1eV_A2_pd.h5')
pdDF, *_ = ep.multiDimXrToPD(data, colDims = 'Eke')

# Pandas .to_hdf - note this also needs a key set, which will be used as the top-level node in the output file.
# Similarly to the above HDF5 case, multiple datasets can be written with different group names (keys).
pdDF.to_hdf(dataFile, key = 'N2')
pdDF.to_hdf(dataFile, key = 'N2v2')

# Pandas to_csv
pdDF.to_csv(dataFile.with_suffix('.csv'))
```

```{code-cell} ipython3
# Read back
import pandas as pd

pdIn = pd.read_hdf(dataFile, key = 'N2')
```

```{code-cell} ipython3
# Convert back to Xarray if desired
import xarray as xr

xrNew = xr.DataArray(pdIn)  # This will produce a flat array
xrRecon, *_ = ep.util.misc.restack(xrNew, refDims = 'matE')   # Restack to specified dataType
xrRecon
```

```{code-cell} ipython3
print(xrRecon.equals(data))  # False...
(xrRecon - data).max()   # But data seems OK
```

```{code-cell} ipython3
# Check file with h5py
hf = h5py.File(dataFile.as_posix(), 'r')
# Top-level groups
print(hf.keys())
hf.close()
```

### Numpy

Raw data can be saved directly as numpy arrays, see [the numpy docs for options and methods](https://numpy.org/doc/stable/user/how-to-io.html?highlight=save#). Note, however, that no metadata will be saved in this case, just the raw arrays.

TODO: wrappers for coords, and metadata/side-car files?

```{code-cell} ipython3
import numpy as np

dataFile = Path(dataPath, 'n2_3sg_0.1-50.1eV_A2.npy')

# For a Numpy ND array from Xarray, use self.to_numpy() or self.data
# Note this is the main array contents only, coords and attrs are discarded
np.save(dataFile, data.to_numpy())
```

```{code-cell} ipython3
# Load and check data
npIn = np.load(dataFile)
print(npIn.dtype, npIn.shape)
(data-npIn).max()
```

## Dictionary methods

For more general/basic IO, Xarray objects can be seralized to native types (this is, in fact, already done on the backend for netCDF and HDF5 datatypes).

Some general routines are provided for this, which wrap Xarray's native `.to_dict()` method with some additional functionality. In particular full serialization of all coordinates is performed (including non-dimensional coords), and some additional wrappers/converters are implemented for problematic data types.

TODO: better wrapping/implementation of `sanitizeAttrsNetCDF()` routine, should be generalised.

```{code-cell} ipython3
# Full dictionary conversion with deconstuctDims
dataDict = ep.util.misc.deconstructDims(data, returnType='dict')

# Deconstruction + split complex data types to Re + Im components
dataDictSplit = ep.util.misc.deconstructDims(data, returnType='dict', splitComplex=True)

# Deconstruction + split complex + clean up attrs for safe types only
dataCleanAttrs, fullAttrs, log = ep.util.xrIO.sanitizeAttrsNetCDF(data)
dataCleanDictSplit = ep.util.misc.deconstructDims(dataCleanAttrs, returnType='dict', splitComplex=True)
```

```{code-cell} ipython3
# Convert back to XR including restacking dimensions
xrRecon = ep.util.misc.reconstructDims(dataDict)
xrRecon
```

```{code-cell} ipython3
:tags: []

# Dictionaries can be handled as usual, e.g. dump to Pickle, JSON
    
import json

# JSON dump to file example
with open(Path(dataPath, 'n2_3sg_0.1-50.1eV_A2_dict.json'), 'w') as f:
    
    # json.dumps(dataDict)  # Fails for complex - not supported by JSON
    
    # OK with split complex to re + im format
    # Note this may still fail in some cases if any unsupported data types are in data.attrs.
    try:
        json.dump(dataDictSplit, f)
        print('JSON OK')
    except:
        print('JSON failed, trying with clean attrs')
        json.dump(dataCleanDictSplit, f)   # Try with cleaned-up attrs
        print('JSON with cleaned attrs OK')

# General JSON complex handling, see notes elsewhere...?
# E.g. https://stackoverflow.com/questions/27909658/json-encoder-and-decoder-for-complex-numpy-arrays
```

```{code-cell} ipython3
# The sanitizeAttrsNetCDF log notes any changes made
log
```

## Pickle

As usual, Pickle can be used as a general method, but the usual caveats apply - it requires data structures and classes to be available on read-in, so can break with library versions, and should not be regarded as archival.

Note, however, that Pickle is the fastest way to save many objects, including complex class objects containing multiple dataarrays. (Other methods to be implemented soon.)

```{code-cell} ipython3
import pickle

# Save a dictionary
with open(Path(dataPath, 'n2_3sg_0.1-50.1eV_A2_dict.pickle'), 'wb') as handle:
    pickle.dump(dataDict, handle, protocol=pickle.HIGHEST_PROTOCOL)
    
    
# Save an Xarray
with open(Path(dataPath, 'n2_3sg_0.1-50.1eV_A2_XR.pickle'), 'wb') as handle:
    pickle.dump(data, handle, protocol=pickle.HIGHEST_PROTOCOL)
    
```

```{code-cell} ipython3
# Read back in to test
with open(Path(dataPath, 'n2_3sg_0.1-50.1eV_A2_dict.pickle'), 'rb') as handle:
    dataDictPklIn = pickle.load(handle)
    
print(dataDict==dataDictPklIn)  # This is False...
data.equals(ep.util.misc.reconstructDims(dataDictPklIn))  # Also False
(ep.util.misc.reconstructDims(dataDictPklIn) - data).max()   # But data looks OK. Dim ordering and/or attrs are different?
```

```{code-cell} ipython3
# Read back in to test
with open(Path(dataPath, 'n2_3sg_0.1-50.1eV_A2_XR.pickle'), 'rb') as handle:
    dataDictPklIn = pickle.load(handle)
    
print(data.equals(dataDictPklIn))  # This is True
(dataDictPklIn - data).max()   # And data looks OK.
```

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
!git ls-remote --heads https://github.com/phockett/ePSproc
```
