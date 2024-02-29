"""
ePSproc conversion functions

28/02/24    Added datasetStack(), mainly for use with multiJob class datastructure.
17/07/20    Added orb3DCoordConv(), for orbital file coord conversions.
12/03/20    Added multiDimXrToPD(), function adapted from lmPlot() code.

TODO: consider implementing CCLIB for unit conversions. See also pint library.

"""

import scipy.constants
import numpy as np
import xarray as xr

from epsproc.util.misc import deconstructDims, reconstructDims, restack
from epsproc.util.selectors import matEleSelector

#************* Xarray handling (to Pandas, to flat, to dict)
def multiDimXrToPD(da, colDims = None, rowDims = None, thres = None, squeeze = True, dropna = True, fillna = False,
                    colRound = 2, verbose = False):
    """
    Convert multidim Xarray to stacked Pandas 2D array, (rowDims, colDims)

    Parameters
    ----------
    da : Xarray
        Array for conversion.

    colDims : list of dims for columns, default = None

    rowDims : list of dims for rows, default = None

    NOTE: either colDims or rowDims must be defined, other coords will be stacked automatically.
    For full control over dim stack ordering, specifiy both colDims and rowDims

    NOTE: if xDim is a MultiIndex, pass as a dictionary mapping, otherwise it may be unstacked during data prep.
    E.g. for plotting stacked (L,M), set xDim = {'LM':['L','M']}

    thres : float, optional, default = None
        Threshold values in output (pd table only)
        TODO: generalise this and use matEleSelector() for input?

    squeeze : bool, optional, default = True
        Drop singleton dimensions.

    dropna : bool, optional, default = True
        Drop all NaN dimensions from output pd data frame (columnwise and rowise).

    fillna : bool, optional, default = False
        Fill any NaN values with 0.0. Useful for plotting/making data contiguous.

    colRound : int, optional, default = True
        Round column values to colRound dp. Only applied for Eke, Ehv, Euler or t dimensions.

    Returns
    -------
    daRestackpd : pandas data frame (2D) with sorted data.

    daRestack : Xarray with restacked data.


    Method
    -------
    Restack Xarray by specified dims, including basic dims checking, then use da.to_pandas().


    12/03/20 Function adapted from lmPlot() code.

    Note
    -----

    This might casue :py:func:`epsproc.lmPlot()` to fail for singleton x-dimensions if squeeze = True. TO do: add work-around, see lines 114-122.

    """

    # Threshold full array - this is as per lmPlot() code, but won't work for xDim as dict.
    # Should set this as a separate function to wrap matEleSelector for general cases.
    # if thres is not None:
    #     # Threshold on abs() value before setting type, otherwise all terms will appear for some cases (e.g. phase plot)
    #     da = matEleSelector(da, thres=thres, inds = selDims, dims = xDim) # , sq = True)  # Squeeze may cause issues here if a singleton dim is used for xDim.


    dimUS = da.unstack().dims

    if type(colDims) == dict:
        colDimsList = list(colDims.items())[0][1]
    elif type(colDims) == list:
        colDimsList = colDims
    else:
        colDimsList = [colDims]  # Workaround for singleton dims not passed as a list - without this set logic below fails.

    if rowDims is None:
        rowDims = list(set(dimUS) - set(colDimsList))   # Use set arithmetic to get items
        rowDims.sort()                      # Set sort to return alphebetical list.


    # Check rowDims exist, otherwise may throw errors with defaults
    rowDimsRed = []
    # for dim in daRestack.unstack().dims:
    #     if dim in rowDims:
    #         rowDimsRed.append(dim)
    for dim in rowDims:
        if dim in dimUS:
            rowDimsRed.append(dim)

    # Additional check for any missing dims
    # Check # of dims and correct for any additional/skipped dims
    # Bit ugly - should be integrated with above code
    if (len(colDimsList) + len(rowDimsRed)) != len(dimUS):
        for dim in dimUS:
            if not (dim in colDimsList) and not (dim in rowDimsRed):
                rowDimsRed.append(dim)

                if verbose:
                    print(f'Adding {dim} to plotting dim list.')


    # Restack for plotting, and drop singleton dimensions if desired.
    # NOTE - plotDim name retained here for compatibility with lmPlot(), may change in future.
    daRestack = da.unstack().stack(plotDim = rowDimsRed).dropna(dim = 'plotDim', how = 'all')

    # Rounding for column values to prevent large float labels in some cases
    for dim in colDimsList:
        if (dim in ['Eke', 'Ehv', 'Euler', 't']) and (colRound is not None):
            daRestack[dim] = daRestack[dim].round(colRound)


    # Restack colDims in cases where it is a MultiIndex
    if type(colDims) == dict:
        daRestack = daRestack.stack(colDims)

    # TODO: add work-around here for singleton x-dim to avoid dropping in that case. (Otherwise have to manually set squeeze = True)
    if squeeze:
        # daRestackpd = daRestack.unstack().stack(plotDim = rowDimsRed).squeeze().to_pandas().dropna(axis = 1).T
        # daRestackpd = daRestack.unstack().stack(plotDim = rowDimsRed).dropna(dim = 'plotDim', how = 'all').squeeze().to_pandas().T
        daRestackpd = daRestack.squeeze().to_pandas()

    else:
        # daRestackpd = daRestack.unstack().stack(plotDim = rowDimsRed).to_pandas().dropna(axis = 1).T
        # daRestackpd = daRestack.unstack().stack(plotDim = rowDimsRed).dropna(dim = 'plotDim', how = 'all').to_pandas().T
        daRestackpd = daRestack.to_pandas()

    # Transpose Pandas table if necessary - colDims must be columns
    if type(colDims) != dict:
        if hasattr(daRestackpd, 'columns') and (colDims not in daRestackpd.columns.names):  # Check coldims, won't exist in singleton dim case.
            daRestackpd = daRestackpd.T

    # For dictionary case, check items for each key are in column names.
    # THIS CODE IS HORRIBLE - should be a neater way to do this.
    # TODO: fix case for some levels missing, at the moment assumes any present is OK.
    # TODO: test for single vs. MultiIndex case - columns.names vs. columns.name?
    else:
        for key in colDims:
            dimList = colDims[key]
            check = [item in daRestackpd.columns.names for item in dimList]
            if not any(check):
                daRestackpd = daRestackpd.T

    # Threshold by abs value, pd only
    # TODO: replace with more general thresholding of input da.
    # AS is this can result in non-contiguous row/col data.
    if thres is not None:
        daRestackpd = daRestackpd[daRestackpd.abs() >= thres]

    # Drop all na cols (often generated by XR restack)
    # NOTE that if data is NOT thresholded, this may do nothing as values may be 0.00, rather than Nan.
    # As set this will drop any all-NaN rows and cols.  Ordering shouldn't matter.
    if dropna:
        daRestackpd = daRestackpd.dropna(how='all', axis=1).dropna(how='all', axis=0)

    # Fill nas for contiguous plotting.
    if fillna:
        daRestackpd = daRestackpd.fillna(0.0)

    return daRestackpd, daRestack


def multiDimXrToDict(da):
    """Convert multiDim Xarray to native dictionary format"""

    # Flatten dims
    daFlat = deconstructDims(da)

    # Convert to dict with native method
    daDict = daFlat.to_dict()

    return daFlat, daDict


def multiDimXrFromDict(daDict):
    """Convert multiDim Xarray to native dictionary format"""

    # Create Xarray
    daFlat = xr.DataArray.from_dict(daDict)

    # Try and rebuild stacked dims either from dataType or dimMap.
    if 'dimMaps' not in daFlat.attrs.keys():

        # Try and restack according to dataType if set
        if 'dataType' in daFlat.attrs.keys():
            da = restack(daFlat)
        else:
            print("*** Can't restack array, missing self.attrs['dimMaps'] and self.attrs['dataType']. For general restacking, try epsproc.util.misc.restack().")
            da = None
    else:
        da = reconstructDims(daFlat)

    return daFlat, daDict


def datasetStack(data, dataType = 'XS', keys = None, stackDim = 'Orb', dimLabel = None,
                 swapDims = False, dropDims = False, sortByStacked = True, **kwargs):
    """
    General routine for stacking dict of XR.dataArray data to new array.

    For use with multiJob class or general data.data[job][dataType] style data. (Cf. TMOdev code...?)

    NOTE: that this format doesn't work for Ehv coords for plotting directly, may need DataSet format for this case.
    NOTE: additonal options for Eke,Ehv dims, which are (orb,E) dependent. May have better routines for handling this elsewhere.

    See also epsproc.util.misc and epsproc.util.xrIO for complementary restacking XR functions.

    Pass **kwargs for matEleSelector to subselect data prior to stacking.


    Parameters
    ----------

    data : dict
        Dictionary of data to stack.
        Should contain elements of type data[key][dataType], e.g. as used by multijob class.

    dataType : str, default = 'XS'
        Data type to stack.

    keys : list, default = None
        List of keys to use.
        If None, skip routine return empty items only.

    stackDim : str, default = 'Orb'
        Dim for concat to new array.
        Can be an existing or new dim.
        Default case assumes data[key] correspond to different orb (channels) jobs.

    dimLabel : str, default = None
        Label for array along stackDim.
        If None, use `data[key][dataType].attrs['jobLabel']` or integer if missing.
        If passed try and use data[key]['jobNotes'][dimLabel].

    swapDims : bool, default = False
        If True swap Eke > Ehv dims.

    dropDims : bool, default = False
        If True, drop secondary Eke or Ehv dim.
        This is currently required for stacking to DA for dims which are (orb,E) dependent.

    sortByStaked : bool, default = True
        If True, run da.sortby() along stackDim.

    **kwargs : optional
        Additional args for ep.matESelector.
        No subselection if not passed.


    Returns
    -------
    xrDA, xrDS, dataDict
        XR dataarray of stacked data.
        XR dataset of stacked data (note this may be None in some cases, when dims are an issue).
        dataDict contains data subset in dictionary format.


    v1  28/02/24  Mainly from quick hacks in https://phockett.github.io/ePSdata/OCS-preliminary/OCS_orbs8-11_AFBLMs_VM-ADMs_140122-JAKE_tidy-replot-200722_v5.html
                  Should check codebase for other complementary routines!

    """

    # Skip empty case.
    if keys is None:
        return None, None, None

    dataDict = {}
    for n,key in enumerate(keys):
        # Set data
        try:
            dataDict[key] = data[key][dataType]
        except KeyError:
            print(f"*** Warning: {key} missing dataType = {dataType}.")
            break

        # Subselect?
        if kwargs:
            dataDict[key] = matEleSelector(data[key][dataType], **kwargs)

        # Test Ehv stacking - may be necessary to maintain this in DA format...?
        # Working with no swap - all orbs have same Eke values.
        # With swap also working, but get many NaNs - issue is that Ehv coords distinct per Orb.
        if swapDims:
            dataDict[key] = dataDict[key].swap_dims({'Eke':'Ehv'})

            if dropDims:
                dataDict[key] = dataDict[key].drop('Eke')

        if dropDims:
            dataDict[key] = dataDict[key].drop('Ehv')


        # Expand dims for stacking
        # TODO: try/except for better logic here for missing labels.
#         dataDict[key].name = dataType
        if stackDim not in dataDict[key].dims:
            if dimLabel is not None:
                if 'jobNotes' in data[key].keys():
                    dataDict[key] = dataDict[key].expand_dims({stackDim:[data[key]['jobNotes'][dimLabel]]})
                else:
                    dataDict[key] = dataDict[key].expand_dims({stackDim:[dataDict[key].attrs[dimLabel]]})

            else:
                if 'jobLabel' in dataDict[key].attrs.keys():
                    dataDict[key] = dataDict[key].expand_dims({stackDim:[dataDict[key].attrs['jobLabel']]})
                else:
                    dataDict[key].expand_dims({stackDim:[n]})


    # Stack to DA... note this may require list not dict (may depend on XR version)
    xrDA = xr.concat([dataDict[k] for k in dataDict.keys()] , stackDim)

    # Sort along new dim?
    if sortByStacked:
        xrDA = xrDA.sortby(stackDim)

    xrDA.name = f"{dataType}"
    xrDA.attrs['jobLabel'] = f"{stackDim}-stacked data"
    xrDA.attrs['stackedKeys'] = keys

    # Stack to DS
    # Should work, but in XR15 get merge issues with Ehv, may just want to reset coord? Should be OK for E-subselected data?
    # MergeError: conflicting values for variable 'Ehv' on objects to be combined. You can skip this check by specifying compat='override'.
    # But this still fails for XR15
    #
    # Running swapDims with dataDict[key].drop('Eke') allows this to work, although format not so useful.
    #
    # TODO: debug and/or check old codes, got this working elsewhere before with additional coord transforms?
    #

    if dropDims:
        xrDS = xr.Dataset(dataDict)
    else:
        xrDS = None

#     xrDA.to_dataset()  # Should also work, but may silently drop some data...? TBC, tested elsewhere already?

#     return xr.concat(dataDict, stackDim)
    return xrDA, xrDS, dataDict



#********************** Calculations

# Convert eV <> Hartrees (atomic units)
def conv_ev_atm(data, to = 'ev'):
    """
    Convert eV <> Hartree (atomic units)

    Parameters
    ----------
    data : int, float, np.array
        Values to convert.

    to : str, default = 'ev'
        - 'ev' to convert H > eV
        - 'H' to convert eV > H

    Returns
    -------
    data converted in converted units.

    """

    # Define H in eV, value from https://en.wikipedia.org/wiki/Hartree_atomic_units#Units
    # H = 27.211386245
    H = scipy.constants.physical_constants['Hartree energy in eV'][0]  # Use scipy value

    if to == 'ev':
        dataOut = data*H
    else:
        dataOut = data/H

    return dataOut

# Convert energy in eV to wavelength in nm
def conv_ev_nm(data): #, to = 'nm'):
    """Convert E(eV) <> lambda(nm)."""

    # Define constants from scipy.constants
    h = scipy.constants.h
    c = scipy.constants.c
    evJ = scipy.constants.physical_constants['electron volt-joule relationship'][0]

    # Define output units - wavelength in m
    waveConv = 1e-9
    dataOut = (h * c)/(data * evJ)/waveConv

    # if to is 'nm':
    #     dataOut = (h * c)/(data * evJ)/waveConv
    #
    # else:
    #     dataOut = (data/waveConv)*evJ/(h * c)

    return dataOut

# Convert energy in eV to wavelength in A, for electrons
def conv_ev_nm_elec(data): #, to = 'nm'):
    """Convert Eke(eV) > lambda(A) for electrons."""

    # Define constants from scipy.constants
    h = scipy.constants.h
    # c = scipy.constants.c
    m = scipy.constants.m_e
    # ec = scipy.constants.e
    evJ = scipy.constants.physical_constants['electron volt-joule relationship'][0]

    # Define output units - wavelength in m
    waveConv = 1e-10
#     dataOut = (h * c)/(data * evJ)/waveConv

    # KE = data*ec
    # nu = np.sqrt(2*KE/m)
    nu = np.sqrt(2*(data*evJ)/m)
    lam = h/(m*nu)/waveConv

    # return (lam, nu, KE)
    return lam

# Renorm by L=0 term
def renormL0(data):
    """
    Renormalise passed data (Xarray) by (L,M) = (0,0) term.

    Requires input Xarray to have dims (L,M) or (l,m), should be robust over all other dims.

    """
    dataOut = data.copy()

    # Note - this currently assumes m dim is present, and forces it to be dropped after selection.
    if hasattr(dataOut,'L'):
        # dataOut /= dataOut.sel({'L':0}).drop('BLM')
        # dataOut /= dataOut.sel({'L':0}).drop('M').squeeze()
        # dataOut = dataOut/dataOut.sel({'L':0}).drop('M').squeeze()  # Non-in-place version, more robust
        dataOut = dataOut/dataOut.sel({'L':0}).sel({'M':0}).drop('M').squeeze()  # Force m=0, issues with spurious m presisting in some cases otherwise

    elif hasattr(dataOut,'l'):
        # dataOut /= dataOut.sel({'l':0}).drop('BLM')
        # dataOut /= dataOut.sel({'l':0}).drop('m').squeeze()
        # dataOut = dataOut/dataOut.sel({'l':0}).drop('m').squeeze()  # Non-in-place version, more robust
        dataOut = dataOut/dataOut.sel({'l':0}).sel({'m':0}).drop('m').squeeze()  # Force m=0, issues with spurious m presisting in some cases otherwise

    else:
        print("***Warning, L/l not present in dataset.")
        return None

    # Propagate attrs
    dataOut.attrs = data.attrs

    return dataOut

# Convert expansion parameters from Legendre Polynomial to Spherical Harmonic form (and back)
def conv_BL_BLM(data, to = 'sph', renorm = True):
    r"""
    Convert BL (Legendre Polynomial) <> BLM (Spherical Harmonic), plus parameter renomalisation.

.. math::
    \beta^{Sph}_{L,0} = \sqrt{(2L+1)/4\pi}\beta^{Lg}

    Note: other conventions may be used here, see https://shtools.github.io/SHTOOLS/complex-spherical-harmonics.html#supported-normalizations

    Parameters
    ----------
    data : Xarray
        Values to convert.
        Currently assumes an Xarray, with dims .L and .M

    to : str, default = 'sph'
        - 'sph' to convert BL > BLM
        - 'lg' to convert BL0 > BL

    renorm : bool, optional, default = True
        If true, additionally renormalise paramters by L=0 term, such that B0 = 1.

    Notes
    -----
    - Should add type to keep track of betas here.
    - Should generalise to other input structure & add error checking.
    - Implement SHTOOLS library....!

    """

    # Set conversion factor
    # Bconv = np.sqrt(2*data.L+1)/(4*np.pi)
    if hasattr(data,'L'):
        Bconv = np.sqrt((2*data.L+1)/(4*np.pi))
    elif hasattr(data,'l'):
        Bconv = np.sqrt((2*data.l+1)/(4*np.pi))

    elif hasattr(data,'XC'):
        # Case for ePS GetCro LF (B2 only) output, with no sigma/XC renorm
        # Q: Sigma renorm in this case...?
        # Bconv = xr.DataArray([np.sqrt(1/(4*np.pi)), np.sqrt(5/(4*np.pi))], dims='XC', coords={'XC':['SIGMA','BETA']})
        Bconv = xr.DataArray([1, np.sqrt(5/(4*np.pi))/np.sqrt(1/(4*np.pi))], dims='XC', coords={'XC':['SIGMA','BETA']})
        renorm = False # No further renorm in this case

    else:
        print("*** Beta conversion error: Data type not supported.")
        return None

    # Set output values
    if to == 'sph':
        dataOut = data/Bconv
    elif to == 'lg':
        dataOut = data*Bconv
    else:
        print(f"*** Beta conversion error: conversion type {to} not supported.")

    if renorm:
        # Note - this currently assumes m dim is present, and forces it to be dropped after selection.
        # if hasattr(dataOut,'L'):
        #     # dataOut /= dataOut.sel({'L':0}).drop('BLM')
        #     dataOut /= dataOut.sel({'L':0}).drop('M').squeeze()
        # elif hasattr(dataOut,'l'):
        #     # dataOut /= dataOut.sel({'l':0}).drop('BLM')
        #     dataOut /= dataOut.sel({'l':0}).drop('m').squeeze()

        # Now moved to separate function
        dataOut = renormL0(dataOut)

    # Propagate attrs
    dataOut.attrs = data.attrs
    dataOut.attrs['normType'] = to

    if 'harmonics' in dataOut.attrs.keys():
        dataOut.attrs['harmonics']['dtype'] = to

    return dataOut

#
def orb3DCoordConv(fileIn, coordMaxLen=50):
    """
    Basic coord parse & conversion for volumetric wavefunction files from ePS.

    Parameters
    ----------
    fileIn : data from a single file
        List of values from a wavefunction file, as returned by :py:func:`epsproc.readOrb3D()`.
        (Note this currently assumes a single file/set of values.)

    coordMaxLen : int, optional, default=50
        Max coord grid size, assumed to demark native Cart (<coordMaxLen) from Spherical (>coordMaxLen) coords.

    Returns
    -------
    x,y,z : np.arrays of Cartesian coords (x,y,z)

    """


    # Set grid, convert to Cart if necessary, assuming that grid won't be larger than 10 Angs
    if (len(fileIn[2][0]) > coordMaxLen):
        # Convert to Cart grid for plotting
        # TODO: Investigate use of sph grid here - should be cleaner.
        # TODO: Investigate recreating mesh in Paraview, rather than saving to file.
        [T,R,P] = np.meshgrid(fileIn[2][1], fileIn[2][0], fileIn[2][2])
        T = (T*np.pi/180) #-np.pi/2
        P = P*np.pi/180
        x = R*np.sin(P)*np.cos(T)
        z = R*np.cos(P)
        y = R*np.sin(P)*np.sin(T)
    else:
        x,y,z = np.meshgrid(file[2][1], file[2][0], file[2][2])

    return x,y,z
