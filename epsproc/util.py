# -*- coding: utf-8 -*-
"""
ePSproc utility functions

Collection of small functions for sorting etc.

11/08/19    Added matEleSelector

"""

import numpy as np

# Selector function for matrix elements in Xarray
def matEleSelector(da, thres = None, inds = None, sq = False, drop=True):
    """
    Select & threshold raw matrix elements in an Xarray

    Parameters
    ----------
    da : Xarray
        Set of matrix elements to sub-select
    thres : float, optional, default None
        Threshold value for abs(matElement), keep only elements > thres.
    inds : dict, optional, default None
        Dicitonary of additional selection criteria, in name:value format.
        These correspond to parameter dimensions in the Xarray structure.
        E.g. inds = {'Type':'L','Cont':'A2'}
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
    if thres is not None:
        daOut = da.where(np.abs(da) > thres, drop = drop)
    else:
        daOut = da

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


# Return list of standard dataArray dims for matrix elements
def matEdimList(sType = 'stacked'):
    """
    Return standard list of dimensions for matrix elements.

    Parameters
    ----------
    sType : string, optional, default = 'stacked'
        Selected 'stacked' or 'unstacked' dimensions.
        Set 'sDict' to return a dictionary of unstacked <> stacked dims mappings for use with xr.stack({dim mapping}).

    Returns
    -------
    list : set of dimension labels.

    """
    if sType is 'stacked':
        # stackedDims
        return ['LM', 'Eke', 'Sym', 'mu', 'it', 'Type']

    elif sType is 'sDict':
        return {'LM':['l','m'],'Sym':['Cont', 'Targ', 'Total']}

    else:
        # unstackedDims
        return ['l','m', 'Eke', 'Cont', 'Targ', 'Total', 'mu', 'it', 'Type']

# Return list of standard dataArray dims for BLM values
def BLMdimList(sType = 'stacked'):
    """
    Return standard list of dimensions for calculated BLM.

    Parameters
    ----------
    sType : string, optional, default = 'stacked'
        Selected 'stacked' or 'unstacked' dimensions.
        Set 'sDict' to return a dictionary of unstacked <> stacked dims mappings for use with xr.stack({dim mapping}).

    Returns
    -------
    list : set of dimension labels.

    """
    if sType is 'stacked':
        # stackedDims
        return ['Euler', 'Eke', 'BLM']

    elif sType is 'sDict':
        return {'BLM':['l','m'],'Euler':['P','T','C']}

    else:
        # unstackedDims
        return ['Eke', 'l', 'm', 'P', 'T', 'C']
