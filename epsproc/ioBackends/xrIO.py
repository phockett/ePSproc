"""
ePSproc wrappers for Xarray file IO.

Native types file IO with Xarray, with additional handling for multiple dimensions and complex numbers.

See the Xarray docs for more backend details: https://docs.xarray.dev/en/stable/user-guide/io.html

See https://github.com/phockett/ePSproc/issues/8 for ongoing notes.

27/06/22    Split out from core IO.py, to extend backend options and support.
            Now additionally wrapped therein for flexible handling of multiple backends.

"""


# Imports
import os
# import re
import numpy as np
# import pandas as pd
# from io import StringIO
import xarray as xr


from epsproc.util.misc import fileListSort, restack
from epsproc.util.xrIO import splitComplex, combineComplex, sanitizeAttrsNetCDF
from epsproc.util.io import setTimeStampedFileName

#**************** Wrappers for Xarray load/save netCDF

# File write wrapper.
def writeXarray(dataIn, fileName = None, filePath = None, engine = 'h5netcdf', forceComplex = False):
    """
    Write file to netCDF format via Xarray .to_netcdf() method.

    Parameters
    -----------
    dataIn : Xarray
        Data array to write to disk.

    fileName : str, optional, default = None
        Filename to use.
        If set to None (default) the file will be written with a datastamp.

    filePath : str, optional, default = None
        Full path to file.
        If set to None (default) the file will be written in the current working directory (as returned by `os.getcwd()`).

    engine : str, optional, default = 'h5netcdf'
        netCDF engine for Xarray to_netcdf method. Some libraries may not support multidim data formats.
        See https://docs.xarray.dev/en/latest/user-guide/io.html

    forceComplex : bool, optional, default = False
        For h5netcdf engine only, set `invalid_netcdf` option = forceComplex.
        If True, complex data will be written directly to file.
        Note this also needs to be read back with the same engine & settings.
        For more details see https://github.com/h5netcdf/h5netcdf#invalid-netcdf-files

    Returns
    -------
    str
        Indicates save type and file path.

    Notes
    -----
    The default option for Xarray is to use Scipy netCDF writer, which does not support complex datatypes. In this case, the data array is written as a dataset with a real and imag component.

    This routine assumes a DataArray as input, although the Xarray file writer pushes this to a DataSet.

    TODO: implement try/except to handle various cases here, and test other netCDF writers (see http://xarray.pydata.org/en/stable/io.html#netcdf).

    Multi-level indexing is also not supported, and must be serialized first. Ugh.

    02/06/22: added improved complex number handling & attibutes sanitizer (lossy).

    """

    if fileName is None:
        # timeString = dt.now()
        # fileName = 'ep_' + timeString.strftime('%Y-%m-%d_%H-%M-%S')
        fileName = setTimeStampedFileName()

    if filePath is None:
        filePath = os.getcwd()

    # Serialize MultiIndex  - testing here for BLM case.
    # if 'BLM' in dataIn.dims:
    #     dataIn = dataIn.reset_index(['Euler','BLM'])

    # Serialize general - use unstact() to flatten all dims
    dataIn = dataIn.unstack()

# Try/except not yet working, multiple error types to handle here...
    # try:
    #     dataIn.to_netcdf(fileName)
    #     saveMsg = 'Written to netCDF4.'
    #     print(saveMsg)
    #     return saveMsg
    #
    # except ValueError as e:
    #     if e.msg != "NetCDF 3 does not support type complex128":
    #         raise
    #     else:
    #         xr.Dataset({'Re':dataIn.real, 'Im':dataIn.imag}).to_netcdf(fileName)
    #         saveMsg = 'Written to netCDF3 (re/im format).'
    #         print(saveMsg)
    #         return saveMsg
    #
    # return 'File not written.'


    if (engine != 'h5netcdf') or (not forceComplex):
        # Safe version with re/im split save type only.
        # Works for scipy and h5netcdf OK, latter will save complex type too, but is strictly not valid.
        dataOut = xr.Dataset({'Re':dataIn.real, 'Im':dataIn.imag})
        # dataOut.attrs = dataIn.attrs   # This will push dataarray attrs to dataset attrs, otherwise they're nested
                                        # May not always want this?

        # Allow for SF & XS coords which may also be complex
        # if 'XS' in dataOut.coords:
        #     dataOut['XSr'] = dataOut.XS.real
        #     dataOut['XSi'] = dataOut.XS.imag
        #     dataOut = dataOut.drop('XS')
        #
        # if 'SF' in dataOut.coords:
        #     dataOut['SFr'] = dataOut.SF.real
        #     dataOut['SFi'] = dataOut.SF.imag
        #     dataOut = dataOut.drop('SF')

        # Allow for arb complex coords.
        # May also want to add attr checker here? Or set in 'sanitizeAttrsNetCDF'
        for item in dataOut.coords.keys():
            if dataOut.coords[item].dtype == 'complex128':
                dataOut.coords[item + 'r'], dataOut.coords[item + 'i'] = splitComplex(dataOut.coords[item])
                dataOut = dataOut.drop(item)

    else:
        # dataOut = dataIn.to_dataset()   # Set to dataset explicitly prior to save - may also need/want to set name here if missing.
        dataOut = dataIn   # Without additional conversion.

    # For netCDF3 can't have multidim attrs, quick fix here for removing them (BLM case)
    if engine == 'scipy':
        if 'sumDims' in dataOut.attrs:
            dataOut.attrs['sumDims'] = [] # test.attrs['selDims'][0]
        if 'selDims' in dataOut.attrs:
            dataOut.attrs['selDims'] = []

    try:
        if engine != 'h5netcdf':
            dataOut.to_netcdf(os.path.join(filePath, fileName + '.nc'), engine=engine)
        else:
            dataOut.to_netcdf(os.path.join(filePath, fileName + '.nc'), engine=engine, invalid_netcdf=forceComplex)

        saveMsg = [f'Written to {engine} format']

    except Exception as e:

        print(f'writeXarray caught exception: {e}')
        print(f'Retrying file write with sanitized attrs.')

        # THIS IS WORKING IN TESTING, but not here?
        # Does work for invalid_netcdf = True case, for h5netcdf backend at least.
        # Seems to be an issue with sanitizeAttrsNetCDF() and/or backend to fix?
        # AH - issue is DataSet vs. DataArray. In former case may need to fix attrs per data variable too.
        # TODO: make use of sanitizeAttrsNetCDF log return. Write to sidecar file?
        dataOut, attrs, log = sanitizeAttrsNetCDF(dataOut)
        if isinstance(dataOut, xr.core.dataset.Dataset):
            for item in dataOut.data_vars:
                dataOut[item], attrs, log = sanitizeAttrsNetCDF(dataOut[item])

        if engine != 'h5netcdf':
            dataOut.to_netcdf(os.path.join(filePath, fileName + '.nc'), engine=engine)
        else:
            dataOut.to_netcdf(os.path.join(filePath, fileName + '.nc'), engine=engine, invalid_netcdf=forceComplex)

        saveMsg = [f'Written to {engine} format, with sanitized attribs (may be lossy)']

    saveMsg.append(os.path.join(filePath, fileName + '.nc'))
    print(saveMsg)

    return saveMsg



# File read wrapper.
def readXarray(fileName, filePath = None, engine = 'h5netcdf', forceComplex = False, forceArray = True):
    """
    Read file from netCDF format via Xarray method.

    Parameters
    -----------
    fileName : str
        File to read.

    filePath : str, optional, default = None
        Full path to file.
        If set to None (default) the file will be written in the current working directory (as returned by `os.getcwd()`).


    engine : str, optional, default = 'h5netcdf'
        netCDF engine for Xarray to_netcdf method. Some libraries may not support multidim data formats.
        See https://docs.xarray.dev/en/latest/user-guide/io.html

    forceComplex : bool, optional, default = False
        For h5netcdf engine only, set `invalid_netcdf` option = forceComplex.
        If True, complex data will be written directly to file.
        Note this also needs to be read back with the same engine & settings.
        For more details see https://github.com/h5netcdf/h5netcdf#invalid-netcdf-files

    forceArray : bool, optional, default = False
        Force file reader to use xr.open_dataarray if True.
        Otherwise use xr.open_dataset.
        The latter case may need additional post-processing, but works with cases with split Re & Im dataarrays.
        If forceComplex = False this setting is ignored.


    Returns
    -------
    Xarray
        Data from file.  May be in serialized format.

    Notes
    -----
    The default option for Xarray is to use Scipy netCDF writer, which does not support complex datatypes. In this case, the data array is written as a dataset with a real and imag component.

    Multi-level indexing is also not supported, and must be serialized first. Ugh.

    TODO: generalize multi-level indexing here.


    07/06/22: improved dim restacking with :py:func:`epsproc.util.misc.restack` routine.
    02/06/22: improved engine & complex number handling (as per writeXarray)

    """

    # TODO - file and path checks
    # See writeXarray() above, and PEMtk.fit._io() functions - should unify method here!
    # See also ep.IO.getFiles?
    # if dataPath is None:
    #     # dataPath = os.getcwd()  # OR
    #     dataPath = Path().absolute()
    # if not Path(fileIn).exists():
    #     fileIn = Path(dataPath,fileIn)  # Assume full path missing if file doesn't exist?

    # Set reader - can try and force to array too.
    # If forceComplex = False, need to use xr.open_dataset for Re+Im dataset format.
    if forceArray and forceComplex:
        freader = xr.open_dataarray
    else:
        freader = xr.open_dataset

    # Read file
    if engine != 'h5netcdf':
        dataIn = freader(fileName, engine = engine)
    else:
        dataIn = freader(fileName, engine = engine, invalid_netcdf = forceComplex)

    if (engine != 'h5netcdf') or (not forceComplex):
        # Reconstruct complex variables, NOTE this drops attrs... there's likely a better way to do this!
        # UPDATE 07/06/22: additional attrs handling below. Note in this case dataOut is a DataArray here.
        dataOut = dataIn.Re + dataIn.Im*1j
        # dataOut.attrs = dataIn.attrs

        # Rest SF & XS coords which may also be complex
        # Note: need to check vs. dataIn here, since dataOut already has dropped vars
        # if 'XSr' in dataIn.data_vars:
        #     dataOut['XS'] = dataIn.XSr + dataIn.XSi*1j
        # #     dataOut = dataOut.drop('XSr').drop('XSi')
        #
        # if 'SFr' in dataIn.data_vars:
        #     dataOut['SF'] = dataIn.SFr + dataIn.SFi
        # #     dataOut = dataOut.drop('SFr').drop('SFi')

        # General version
        for item in dataOut.coords.keys():
            # Check for r+i pairs - note labelling assumed to match writeXarray conventions here.
            if item.endswith('r'):
                itemi = item[:-1] + 'i'

                # If imag partner found, restack and remove split components.
                if itemi in dataOut.coords.keys():
                    dataOut.coords[item[:-1]] = combineComplex(dataOut.coords[item], dataOut.coords[itemi])
                    dataOut = dataOut.drop([item,itemi])

    else:
        dataOut = dataIn

    # For dataset case, try some generic handling. May need more sophisticated methods here, maybe just assume DataArray and convert?
    if (not dataOut.attrs) and isinstance(dataIn, xr.core.dataset.Dataset):
        dataOut.attrs = dataIn[list(dataIn.data_vars)[0]].attrs

    # Recreate MultiIndex from serialized version  - testing here for BLM case.
    # if 'BLM' in dataIn.dims:
    #     dataIn = dataIn.set_index({'BLM':['l','m'],'Euler':['P','T','C']})

    # Recreate MultiIndex from serialized version according to array type.
    # 01/06/22: added try/except for lazy dim handling.
    try:
        # if dataIn.dataType == 'BLM':
        #     dataOut = dataOut.stack(BLMdimList(sType = 'sDict'))
        # elif dataIn.dataType == 'matE':
        #     dataOut = dataOut.stack(matEdimList(sType = 'sDict'))

        dataOut, dims = restack(dataOut)  # General restacking routine, may want to pass args here for more flexibility.

    except:
        print(f"Failed to restack input dataset for dataType {dataIn.dataType}, dims may be missing. Check ep.dataTypesList['{dataIn.dataType}'] for details.")


    return dataOut
