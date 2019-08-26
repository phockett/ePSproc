# -*- coding: utf-8 -*-
"""
ePSproc utility functions

Collection of small functions for sorting etc.

11/08/19    Added matEleSelector

"""

import numpy as np

# Selector function for matrix elements in Xarray
def matEleSelector(da, thres = None, inds = None):
    """
    Select & threshold raw matrix elements in an Xarray

    Parameters
    ----------
    da : Xarray
        Set of matrix elements to sub-select
    thres : float, optional, default None
        Threshold value for abs(matElement), keep only elements > thres.
    inds : dict, optional
        Dicitonary of additional selection criteria, in name:value format.
        These correspond to parameter dimensions in the Xarray structure.
        E.g. inds = {'Type':'L','Cont':'A2'}

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


    # Reduce dims by thesholding on abs values
    # Do this after selection to ensure Nans removed.
    if thres is not None:
        daOut = da.where(np.abs(da) > thres, drop=True)
    else:
        daOut = da

    return daOut


# Select over vals from data structure
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
