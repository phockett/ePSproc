"""
Density matrix routines

27/08/21    v2 Updated dim handling for renaming multi-index levels + working HV plotting routines.
26/08/21    v1 Initial implementation

Dev code:

- http://100.86.127.24/jupyter/user/paul/doc/tree/github/ePSproc/notebooks/in_progress/ePSdev-PEMtk_correlations_den-mats_fn-def_260821.ipynb
- http://100.86.127.24/jupyter/user/paul/doc/tree/github/ePSproc/notebooks/in_progress/ePSdev-PEMtk_correlations_den-mats_basic-tests_220821.ipynb

"""

#*** Dim functionality (see also lmPlot() and multiDimXrToPD() functions)

# Set imports
from epsproc.util import matEleSelector
from epsproc.util.misc import checkDims
# matEdimList, BLMdimList, dataTypesList, multiDimXrToPD
# checkDims = ep.util.misc.checkDims

def dimRestack(da, stackDims = []):
    """
    General Xarray restacker including multi-indexes.

    Check dims in da, and restack according to refDims if necessary

    Parameters
    ----------
    da : xarray
        Data array to check & restack

    stackDims : str, list, dict, optional, default = []
        Dimensions to check in da, and restack along if not already stacked.
        Note that this can mix stacked and unstacked dims, and will restack if necessary

    Returns
    -------
    daOut : Xarray
        Data array with restacked dim.

    stackedDim : str
        Name of new stacked dim

    rsMap : dict
        Dictionary mapping for new stacked dim

    dimCheck : dict
        Full output from :py:func:`ep.util.misc.checkDims()`


    Examples
    --------
    >>> # Assuming matE is a standard array of matrix elements
    >>> daOut, stackedDim, rsMap, dimCheck = dimRestack(matE)  # OK, returns input + dim check results
    >>> daOut, stackedDim, rsMap, dimCheck = dimRestack(matE, stackDims='LM')  # OK, returns original dims
    >>> daOut, stackedDim, rsMap, dimCheck = dimRestack(matE, stackDims=['LM','mu'])  # OK, restacks dims


    Notes
    -----

    Passing new mappings as stackDims is currently not supported, e.g. dMap = {'NewDim':[d1,d2...]} will fail.
    For this case, just use the native da.stack(dMap).


    See also
    --------

    multiDimXrToPD() : map input array to 2D Pandas DataFrame

    """


    #*** Check dims
    dimCheck = checkDims(da, refDims = stackDims)

    if dimCheck['missing']:
        print(f"***Error: Missing specified dimension(s): {dimCheck['missing']}")

        return dimsIn

    #*** Restack dims if necessary
    # (1) if there are both stacked and unstacked dims in denDims, unstack then restack
    if dimCheck['stackedShared'] and dimCheck['shared']:

        # For restack, check stacked vs. unstacked dims and assign new mapping
        rsDims = [dimCheck['stackedMap'][key] for key in dimCheck['stackedShared']]
        rsDims.append(dimCheck['shared'])  # Append unstacked dims

        rsDims = [item for sublist in rsDims for item in sublist]  # Force to 1D list, https://stackoverflow.com/a/952952

        # Set new map
        rsMap = {','.join(rsDims):rsDims}  # Assumes dim names are str, otherwise use ','.join(map(str, rsDims))

        # Restack
        daOut = da.unstack(dimCheck['stackedShared']).stack(rsMap)
        stackedDim = next(iter(rsMap))   # list(rsMap.keys())[0]  # Set name for later

    # (2) for multiple single dims, restack only
    elif len(dimCheck['shared'])>1:

        rsDims = dimCheck['shared']

        # Set new map
        rsMap = {','.join(rsDims):rsDims}  # Assumes dim names are str, otherwise use ','.join(map(str, rsDims))

        # Restack
        daOut = da.stack(rsMap)
        stackedDim = next(iter(rsMap))   # list(rsMap.keys())[0]  # Set name for later

    # (3) if the dim is already stacked, just pass back.
    else:
        rsMap = {}
        daOut = da
        stackedDim = stackDims  # Use passed dim name

    return daOut, stackedDim, rsMap, dimCheck





def densityCalc(da, denDims = 'LM',
                sumDims = None, keepDims = None,   # May want additional control here, 1/sumDims ?
                selDims = None, thres = None, squeeze = True):

    r"""
    General density matrix from Xarray.

    Compute density matrix as (outer product) da[denDims]*da[denDims].conj(), where dim specifies the dimension(s) to use.
    This is, essentially, the density matrix :math:\row = |denDims\rangle \langle denDims|:math:

    Parameters
    ----------
    da : xarray
        Data array to check & restack

    denDims : str, list, dict, optional, default = 'LM'
        Dimensions to use as "state vector" from da.
        If a single dim (including stacked dims), which exists, this will be used directly.
        If multiple dims, will be restacked to a new dimension. Note that this can mix stacked and unstacked dims, and will restack if necessary.

    sumDims : str, list, bool, optional, default = None
        Set specific dims to sum ("trace") over.
        If sumDims = True all dims, apart from denDims and keepDims, will be summed over.

    keepDims : str, list, optional, default = None
        Define dims to keep (won't be summed over). Only used if sumDims = True.

    selDims : str, list, optional, default = None
        Dimensions to subselect from.

    thres : float, optional, default = None
        Threshold value. If set, used for both input and output datasets.

    squeeze : bool, optional, default = True
        Squeeze out singleton dims if True.

    Notes
    -----
    selDims, thres and squeeze are passed to the standard :py:func:`matEleSelector` function.

    """

    # Set data
    daDen = matEleSelector(da, thres = thres, inds = selDims, sq = squeeze)  #.sum(sumDims)
                                                # TODO: pass **kwargs here?
                                                # Pass dims = denDims?

    # Restack dims if required
    daDen, denDim, rsMap, dimCheck = dimRestack(daDen, stackDims = denDims)


    # Check for summation dims, cf. padPlot() routine
#     extraDims = set(daDen.dims) - {*facetDimsCheck,*sumDims}  # Check for outstanding dims, this will return an empty set if all dims accounted for here

#         if extraDims:
#             print(f"Found additional dims {extraDims}, summing to reduce for plot. Pass selDims to avoid.")

#             for dim in extraDims:
#                 subset = subset.sum(dim)  #.squeeze()

    # Handle dim summation from input args only.
    # Default is to leave all other dims untouched
    if sumDims or keepDims:

        # Set this to sum over all other dims (except keepDims)
        if sumDims is True:
            sumDims = {*daDen.dims} - {denDim}

        elif not isinstance(sumDims, list):  # If passed, force list
            sumDims = [sumDims]

        if not isinstance(keepDims, list):  # If passed, force list
            keepDims = [keepDims]

        sumDimsCheck = set(daDen.dims)&{*sumDims} # This checks sumDims are present, otherwise will throw an error.
        keepDimsCheck = set(daDen.dims)&{*keepDims}

        sumTot = sumDimsCheck - keepDimsCheck  # Set dims to sum
        print(sumDimsCheck, keepDimsCheck, sumTot)

    else:
        sumTot = []


    #*** Compute density from specified dims
#     daDot = daDen * daDen.conj().rename({denDim:denDim+'_p'})  # Outer product along denDims
                                                                # This works, but can give issues with shared multi-index dims


    # Version with renaming of multi-index dims prior to outer-product - avoids linked dims in output array.
    newDims = {item:item+'_p' for item in rsMap[denDim]}
    daConj = daDen.conj().unstack(denDim).rename(newDims).stack({denDim+'_p':list(newDims.values())})
    daDot = daDen * daConj

    # General .dot version... might be faster than above?
#     matEdot = xr.dot(daDen, daDen.conj().rename({denDim:denDim+'_p'}), dims = sumDims)

    daOut = matEleSelector(daDot, thres = thres, sq = squeeze).sum(sumTot) # Threshold density mat
                                                                        # TODO: pass **kwargs here?
                                                                        # Pass dims = denDims?

    return daOut, daDot

import pandas as pd
import holoviews as hv
# import hvplot.pandas
hv.extension('bokeh')

def matPlot(da, kdims = ['LMp','LM'], pTypes = ['a','i','r'], returnType = 'plot'):
    """
    General matrix (2D) plot + stacked dims plotter with HoloViews.

    """

    #*** Set data
    daPlot = da.copy()  # May want to add thresholds etc. here?
    attrs = daPlot.attrs.copy()

    #*** Relabel any multi-indexes to strings.
    # Without this HV plot throws errors for int category dims (kdims only?)
    # TODO: check for kdims vs. stacked/selection dims?
    ind = daPlot.indexes

    labels = {}
    for dim in ind:
    #     print(ind[item])

        # Remap multiindex only - note this returns list not PD indexer
        if isinstance(ind[dim], pd.core.indexes.multi.MultiIndex):
            labels[dim] = [','.join(map(str, item)) for item in ind[dim]]

            daPlot = daPlot.assign_coords({dim:labels[dim]})

    # Basic version - working for any dim dataset, but single pType only
#     daPlot = ep.plotTypeSelector(matEplot, pType = pType)
#     hvds = hv.Dataset(daPlot)
#     hvmap = hvds.to(hv.HeatMap, kdims=kdims)

    # Compose HoloMap from data for r,i,a plottypes
    # This only works for 2D data...? (Or case where dims are otherwise specified in loop.)
#     hvmap = hv.HoloMap({pType: hv.Dataset(ep.plotTypeSelector(matEplot, pType = pType).real).to(hv.HeatMap, kdims=kdims)
#                          for pType in ['r','i','a']}, kdims="pType")

    # Example with continuum included too - should be able to generalise?
#     hvmap = hv.HoloMap({(cont,pType): hv.Dataset(ep.plotTypeSelector(matEplot.sel(Cont=cont), pType = pType).squeeze(drop=True).real).to(hv.HeatMap, kdims=['LMp','LM'])
#                          for pType in ['r','i','a'] for cont in ['A2','B1','B2']}, kdims=['Cont',"pType"])

    # Stack to xr.Dataset for pTypes...
    # NOTE .copy() here, otherwise end up with null valued output (overwrites/sums?)
    daPlotDS = xr.Dataset({pType:ep.plotTypeSelector(daPlot.copy(), pType = pType) for pType in pTypes}) #['r','i','a']})
    daPlot = daPlotDS.to_array().rename({'variable':'pType'})  # Restack pType to array
    daPlot.attrs = attrs  # Propagate attrs
    daPlot.name = da.name
    hvds = hv.Dataset(daPlot)
    hvmap = hvds.to(hv.HeatMap, kdims=kdims)

    if returnType is 'full':
        return hvmap, hvds, daPlotDS
    else:
        return hvmap
