"""
ePSproc conversion functions

12/03/20    Added multiDimXrToPD(), function adapted from lmPlot() code.

"""

import scipy.constants

def multiDimXrToPD(da, colDims = None, rowDims = None, thres = None, squeeze = True, dropna = True, verbose = False):
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
        Drop all NaN dimensions from output pd data frame (columnwise).

    Returns
    -------
    daRestackpd : pandas data frame (2D) with sorted data.

    daRestack : Xarray with restacked data.


    Method
    -------
    Restack Xarray by specified dims, including basic dims checking, then use da.to_pandas().


    12/03/20 Function adapted from lmPlot() code.

    """

    # Threshold full array - this is as per lmPlot() code, but won't work for xDim as dict.
    # Should set this as a separate function to wrap matEleSelector for general cases.
    # if thres is not None:
    #     # Threshold on abs() value before setting type, otherwise all terms will appear for some cases (e.g. phase plot)
    #     da = matEleSelector(da, thres=thres, inds = selDims, dims = xDim) # , sq = True)  # Squeeze may cause issues here if a singleton dim is used for xDim.


    dimUS = da.unstack().dims

    if type(colDims) == dict:
        colDimsList = list(colDims.items())[0][1]
    else:
        colDimsList = colDims

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

    # Restack colDims in cases where it is a MultiIndex
    if type(colDims) == dict:
        daRestack = daRestack.stack(colDims)


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
        if colDims not in daRestackpd.columns.names:
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
    if thres is not None:
        daRestackpd = daRestackpd[daRestackpd.abs() >= thres]

    # Drop all na cols (often generated by XR restack)
    if dropna:
        daRestackpd = daRestackpd.dropna(how='all', axis=1)


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
    """Convert E(eV) <> nu(nm)."""

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
