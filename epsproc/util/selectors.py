#*************** Selection functions

import numpy as np

# Selector function for matrix elements in Xarray
def matEleSelector(da, thres = None, inds = None, dims = None, sq = False, drop=True):
    """
    Select & threshold raw matrix elements in an Xarray

    Parameters
    ----------
    da : Xarray
        Set of matrix elements to sub-select
    thres : float, optional, default None
        Threshold value for abs(matElement), keep only elements > thres.
        This is *element-wise*.
    inds : dict, optional, default None
        Dicitonary of additional selection criteria, in name:value format.
        These correspond to parameter dimensions in the Xarray structure.
        E.g. inds = {'Type':'L','Cont':'A2'}
    dims : str or list of strs, dimensions to look for max & threshold, default None
        Set for *dimension-wise* thresholding. If set, this is used *instead* of element-wise thresholding.
        List of dimensions, which will be checked vs. threshold for max value, according to abs(dim.max) > threshold
        This allows for consistent selection of continuous parameters over a dimension, by a threshold.
    sq : bool, optional, default False
        Squeeze output singleton dimensions.
    drop : bool, optional, default True
        Passed to da.where() for thresholding, drop coord labels for values below threshold.

    Returns
    -------
    daOut
        Xarray structure of selected matrix elements.
        Note that Nans are dropped if possible.

    Example
    -------
    >>> daOut = matEleSelector(da, inds = {'Type':'L','Cont':'A2'})

    """

    # Iterate over other selection criteria
    # This may return view or copy - TBC - but seems to work as expected.
    # http://xarray.pydata.org/en/v0.12.3/indexing.html#copies-vs-views
    if inds is not None:
        da = da.sel(inds)    # Fors inds as dict, e.g. {'Type':'L','it':1,'Cont':'A2'}
                                # May want to send as list, or automate vs. dim names?
                                # NOTE - in current dev code this is used to reindex, so .squeeze() casuses issues!


    # Reduce dims by thesholding on abs values
    # Do this after selection to ensure Nans removed.
    if (thres is not None) and (dims is None):
        daOut = da.where(np.abs(da) > thres, drop = drop)
    else:
        daOut = da

    # If dims is set, check over dims for consistency.
    # WILL this just produce same results as thres then squeeze...?
    if (dims is not None) and (thres is not None):
        daOut = daOut.where(np.abs(da).max(dim = dims) > thres, drop = drop)

    if sq:
        daOut = daOut.squeeze()  # Squeeze dims.

    return daOut


# Select over vals from data structure (list)
# Currently only used in IO.matEleGroupDim
def dataGroupSel(data, dInd):
    a = data[0]
    dataSub = []

    uVals = np.unique(a[dInd,:])

    for val in uVals:
        # Get matching terms and subset data
        # iSel = np.nonzero(a[dInd,:]==val)
        iSel = (a[dInd,:]==val)
        dataSub.append([data[0][:,iSel], data[1][iSel]])

    return dataSub


# Xarray groupby + compare values
# STARTED... but not finished.  For basic diff along a dimension, just use da.diff(dim), see http://xarray.pydata.org/en/stable/generated/xarray.DataArray.diff.html#xarray.DataArray.diff
# def groupCmp(data, dim):
#     """
#     Basic routine to compare sets of values by dimension, using Xarray groupby functionality.
#
#     Parameters
#     ----------
#     data : Xarray
#         Data for comparison
#
#     dim : str
#         Dimension label for grouping
#
#     Returns
#     -------
#
#     """
#
#     dGroup = data.groupby(dim)
#
#     # Check differences between groups
#     for gTest in dGroup:
