"""
ePSproc setmatE functions

Manual setting for matrix elements.

25/04/23    v1, based on existing setBLMs and setADMs functions.
            Dev work: code-share/jupyter-shared/ePSproc_dev_2022/util_funcs/setMatE_dev_230423.ipynb


TODO:

- Implment conformDims directly in ePSproc, see PEMtk.toePSproc() for method.
- General handling for additional dims, currently variable from data, but only 'Eke' supported as additional mapping. Should change to dict mapping.


"""


import numpy as np
import xarray as xr
import pandas as pd

try:
    from pemtk.sym._util import toePSproc
    toePSprocFlag = True
except ImportError:
    print("*** Missing PEMtk, can't import `toePSproc` method.")
    toePSprocFlag = False

from epsproc.util import listFuncs, misc
from epsproc.util.conversion import multiDimXrToPD
from epsproc.sphFuncs.sphConv import checkSphDims

# General BLM setter for using with custom values.
# 04/04/23: hacking in as per existing setADMs() and cf. also blmXarray().
# TODO: implement dim remapping, see PEMtk.toePSproc()
# v2: moding for variable coords
def setMatE(data = [0,0,0,1], dataNames = ['l','m','mu'], LMLabels = None, symLabels = None, typeLabels = None,   # keyDims = {},
            Evals = None, eType = 'Eke', eUnits = 'eV',
            name = None, conformDims = False, **kwargs):
    """
    Create Xarray from BLMs, or create default case BLM = [0,0,1].

    Parameters
    ----------
    data : list or np.array, default = [0,0,0,1]
        Set of matrix elements = [l,m,mu, value].
        If multiple values are provided per (l,m,mu) index, they are set to the E or t axis (if provided), or indexed numerically.
        For more control, specify either dataName = [l,m,...] to specify dimensions, or pass unlabeled data and set LMLabels list.
        E.g. setmatE(data = [0,0,1], dataNames = ['l','m']) will set only (l,m) values.
             setmatE(data = [1], dataNames = ['l','m'], LMLabels=[0,0]) will set as per the default case with different input list styles.

    dataNames : list, default = ['l','m','mu']
        Dimension names corresponding to input data.

    LMLabels : list or np.array, optional, default = None
        If passed, assume data elements are unabelled, and use indicies provided here (with names from dataNames).

    symLabels : not implemented

    typeLabels : not implemented

    Evals : list or np.array, optional, default = None
        If passed, use for dimension defining data sets (usually time or energy).
        Defaults to numerical label if not passed, Evals = np.arange(0,data.shape[1])

    eType : str, optional, default = 'Eke'
        Energy type, usually Eke or Ehv

    eUnits : str, optional, default = 'eV'
        Energy units.

    name : str, optional, default = None
        Set a name for the array.
        If None, will be set to 'ADM' (same as dataType attrib)

    conformDims : bool, optional, default = False
        Add any missing dims to match ep.listFuncs.BLMdimList.
        Currently requires PEMtk.toePSproc() for function.


    Returns
    -------
    matE : Xarray
        matE in Xarray format, dims as per :py:func:`epsproc.utils.BLMdimList()`

    TODO:

    - Implment conformDims directly in ePSproc, see PEMtk.toePSproc() for method.
    - General handling for additional dims, currently variable from data, but only 'Eke' supported as additional mapping. Should change to dict mapping.

    Examples
    ---------
    >>> # Default case
    >>> matE = setMatE()
    >>> matE

    >>> # Set with custom dims
    >>> matE = setmatE(data=[0,0,1], dataNames=['l','m'])
    >>> matE

    >>> # Set with ranges (as an array [dim1,dim2,..., E(0), E(1)...]), with values per E
    >>> EPoints = 10
    >>> matE = setmatE(data = [[0,0, *np.ones(10)], [2,0, *np.linspace(0,1,EPoints)], [4,0, *np.linspace(0,0.5,EPoints)]], dataNames=['l','m'])
    >>> matE.attrs['PD']   # Show Pandas table

    """

    # Check size of passed set of ADMs
    # For ease of manipulation, just change to np.array if necessary!
    if isinstance(data, list):
        data = np.array(data, ndmin = 2)

    # Force additional dim for 1D case
    elif isinstance(data, np.ndarray):
        if data.ndim < 2:
            data = data[:,np.newaxis]

    # Set lables explicitly if not passed, and resize ADMs
    if LMLabels is None:
        LMLabels = data[:,0:len(dataNames)]
        data = data[:,len(dataNames):]


    # Set indexing, default to numerical
    if Evals is None:
        Evals = np.arange(0,data.shape[1])
        eUnits = 'Index'

    # Set up Xarray
    # return LMLabels

    QNs = pd.MultiIndex.from_arrays(LMLabels.real.T.astype('int8'), names = dataNames)  # Set lables, enforce type

    # return QNs, LMmu
    stackDim = ''.join(dataNames)
    matE = xr.DataArray(data, coords={stackDim:QNs,eType:Evals}, dims = [stackDim,eType])

    # Metadata
    if name is None:
        name = 'matE'
    # else:
    #     matE.name = name

    matE.name = name
    matE.attrs['dataType'] = 'matE'
    matE.attrs['long_name'] = 'Matrix elements (manually defined)'

    attrsDict = listFuncs.dataTypesList()['matE']  # Get standard defn.
    # getRefDims(data = None, refType = None, sType = 'sDict')
    matE.attrs['source'] = 'epsproc.setMatE'
    matE.attrs['desc'] = 'Raw photoionization matrix elements.'
    matE.attrs['def'] = attrsDict['def']


# listFuncs.dataTypesList() currently defines as below... may want to add manual version?
#     {'source': 'epsproc.IO.readMatEle',
#  'desc': 'Raw photoionization matrix elements from ePS, DumpIdy command and file segments.',
#  'recordType': 'DumpIdy',
#  'dims': {'LM': ['l', 'm'], 'Sym': ['Cont', 'Targ', 'Total']},
#  'def': <function epsproc.util.listFuncs.matEdimList(sType='stacked')>}


#     # Set units
#     BLMX.attrs['units'] = 'arb'
#     BLMX.t.attrs['units'] = tUnits

    # Set defaults for harmonics
    matE.attrs['harmonics'] = listFuncs.YLMtype(**kwargs)

    if conformDims:
        if toePSprocFlag:
            # print("*** ConformDims NOT YET IMPLEMENTED - see PEMtk.toePSproc() for method.")
            matE = toePSproc(matE.to_dataset(name='XR'), dataType = 'matE', dimMap={})   # Note - this currently expects a dataset, should modify.
            matE.name = name  # Rename output
            # matE.attrs['dims'] = attrsDict['dims']  # Set to standard mapping?

        else:
            print("*** Missing PEMtk library, can't import `toePSproc` method for `conformDims` functionality.")

        checkSphDims(matE)  # Set sph dim info

        # return matEout

    # else:
    #     return matE

    else:
        checkSphDims(matE, keyDims={stackDim:dataNames})   # Set sph dim info for custom dim mappings

    matE.attrs['dims'] = misc.checkDims(matE)  # Dim mappings info

    # Set PD table
    matE.attrs['PD'],_ = pdTest, _ = multiDimXrToPD(matE, colDims = eType if len(Evals)>1 else matE.attrs['harmonics']['keyDims'], thres=None, squeeze=False)

    return matE
