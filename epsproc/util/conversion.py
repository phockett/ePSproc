"""
ePSproc conversion functions

17/07/20    Added orb3DCoordConv(), for orbital file coord conversions.
12/03/20    Added multiDimXrToPD(), function adapted from lmPlot() code.

TODO: consider implementing CCLIB for unit conversions. See also pint library.

"""

import scipy.constants
import numpy as np
import xarray as xr

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

    if to is 'ev':
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
        dataOut /= dataOut.sel({'L':0}).drop('M').squeeze()
    elif hasattr(dataOut,'l'):
        # dataOut /= dataOut.sel({'l':0}).drop('BLM')
        dataOut /= dataOut.sel({'l':0}).drop('m').squeeze()
    else:
        print("***Warning, L/l not present in dataset.")
        return None

    # Propagate attrs
    dataOut.attrs = data.attrs

    return dataOut

# Convert expansion parameters from Legendre Polynomial to Spherical Harmonic form (and back)
def conv_BL_BLM(data, to = 'sph', renorm = True):
    """
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
    if to is 'sph':
        dataOut = data/Bconv
    elif to is 'lg':
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
