# TODO: HDF5 wrappers for dict forms.
#       Integration with IO.py file IO.
#       See also PEMtk.fit._io for similar functionality/prototypes.
#       And dev in http://jake/jupyter/user/paul/doc/tree/code-share/jupyter-shared/PEMtk_dev_2022/io_dev/xarray_dict_IO_tests_210622-Jake.ipynb


# Imports
import os
# import re
# import numpy as np
# import pandas as pd
# from io import StringIO
import xarray as xr
from pathlib import Path

from ast import literal_eval

try:
    import h5py
    h5Flag = True

except ImportError as e:
    if e.msg != "No module named 'h5py'":
        raise
    print('* h5py not found, HDF5 export not available. ')
    h5Flag = False


from epsproc.util.misc import deconstructDims, reconstructDims
from epsproc.util.xrIO import splitComplex, combineComplex, sanitizeAttrsNetCDF
from epsproc.util.io import setTimeStampedFileName


def writeXarrayToHDF5(data, fileName = None, filePath = None, dataName = None, appendFlag = True):
    """Write Xarray or dictionary to HDF5"""

    if not h5Flag:
        print('*** Install h5py for HDF5 export. ')
        return None

    if filePath is None:
        # dataPath = os.getcwd()  # OR
        filePath = Path().absolute()

    # Create file if None
    if fileName is None:
        fileName = setTimeStampedFileName(ext='h5')

    # Name group
    # TODO: try and get this from input data structure
    if dataName is None:
        dataName = 'xr'

    # Open file, append if existing
    fOut = Path(filePath,fileName)

    # This needs more thought - fails if file exists and data key exists in both cases.
    if fOut.is_file() and appendFlag:
        hf = h5py.File(fOut.as_posix(), 'a')
    else:
        hf = h5py.File(fOut.as_posix(), 'w')

    dict_group = hf.create_group(dataName)

    # Convert Xarray to safe dictionary format
    # TODO: test with DataSets
    if isinstance(data, xr.core.dataarray.DataArray):   # or isinstance(data, xr.core.dataarray.DataSet):
        dataDict = deconstructDims(data, returnType = 'dict')

        # Reformat output for HDF5
        # See also
        for k,v in dataDict.items():
            # print(k)
            # v = dict_test[k]
            # if not v:
            #     dict_group[k] = str(v)
            # else:
            #     dict_group[k] = v

            # This is OK for general handling
            # try:
            #     dict_group[k] = v
            # except:
            #     dict_group[k] = str(v)

            # Force all groups to str except data
            # This avoids issues with tuples, empty items and nested dicts
            # See also ep.util.xrIO.sanitizeAttrsNetCDF() - may want to implement that instead?
            if k == 'data':
                dict_group[k] = v

            else:
                dict_group[k] = str(v)

    hf.close()

    return fOut


def readXarrayFromHDF5(fileName, filePath = None, dataName = None, evalStrings = True):
    """
    Read Xarray or dictionary to HDF5

    TODO: better handling of coords & attrs. Currently only works with string eval, which is hacky.

    """

    if not h5Flag:
        print('*** Install h5py for HDF5 read. ')
        return None

    if filePath is None:
        # dataPath = os.getcwd()  # OR
        filePath = Path().absolute()

    # Name group
    # TODO: try and get this from input data structure
    # TODO: default to reading all?
    if dataName is None:
        dataName = 'xr'

    # Open file
    fIn = Path(filePath,fileName)
    hf = h5py.File(fIn.as_posix(), 'r')

    # Load data - currently assumes single object
    dict_new = {}
    dict_group_load = hf[dataName]
    dict_group_keys = dict_group_load.keys()

    for k in dict_group_keys:
        v= dict_group_load[k][()]   #[:]   # Get data (not just object), see https://docs.h5py.org/en/stable/high/dataset.html#reading-writing-data

    #     print(v)

        # Try converting items if necessary
        # ok - WORKS FOR ALL CASES EXCEPT NON-EXECUTABLE STRS
        # Note this assumes items are safe to eval!
        try:

            if isinstance(v,bytes):
                dict_new[k] = v.decode("utf-8")
            else:
                dict_new[k] = v

            if evalStrings:
                dict_new[k] = literal_eval(dict_new[k])

        # Push to output directly
        except:
            dict_new[k] = v

    hf.close()

    # Rebuild
    # TODO: add some error checking here
    # TODO: also make this optional and add restack() routine?
    try:
        xrFromDict = reconstructDims(dict_new)

        return (dict_new, xrFromDict)

    except:
        print(f'*** Failed to rebuild Xarray from {fIn}, returning dict only.')
        return (dict_new)
