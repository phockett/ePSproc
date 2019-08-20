# -*- coding: utf-8 -*-
"""
ePSproc utility functions

Collection of small functions for sorting etc.

11/08/19    Added matEleSelector

"""


# Selector function for matrix elements in Xarray
def matEleSelector(da, thres = 1e-2, inds = None):
    """
    Select & threshold raw matrix elements in an Xarray
    
    Parameters
    ----------
    da : Xarray
        Set of matrix elements to sub-select
    thres : float, optional, default = 1e-2
        Threshold value for abs(matElement), keep only elements > thres.
    inds : dict, optional
        Dicitonary of additional selection criteria, in name:value format.
        These correspond to parameter dimensions in the Xarray structure.
        E.g. inds = {'ip':1,'Cont':'A2'}
    
    Returns
    -------
    daOut
        Xarray structure of selected matrix elements.
        Note that Nans are dropped if possible.
        
    Example
    -------
    >>> daOut = matEleSelector(da, inds = {'ip':1,'Cont':'A2'})

    """
    
    # Iterate over other selection criteria
    # This may return view or copy - TBC - but seems to work as expected.
    # http://xarray.pydata.org/en/v0.12.3/indexing.html#copies-vs-views
    if inds is not None:
        da = da.sel(inds)    # Fors inds as dict, e.g. {'ip':1,'it':1,'Cont':'A2'}
                                # May want to send as list, or automate vs. dim names?

                                
    # Reduce dims by thesholding on abs values
    # Do this after selection to ensure Nans removed.
    daOut = da.where(np.abs(da) > thres, drop=True)
                                    
    return daOut
    
