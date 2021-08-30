"""
Density matrix routines

30/08/21    v3 Debugged and updated docstrings.
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
                selDims = None, thres = None, squeeze = True,
                keepNaNs = False):

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

    keepNaNs : bool, optional, default = False
        Keep NaNs in result if false.
        Otherwise, set to zero. (Note this is currently implemented indirectly via xr.sum(skipna = ~keepNaNs).)
        Note: in some cases this may cause NaNs to propagate and give all-NaN results.


    Notes
    -----
    selDims, thres and squeeze are passed to the standard :py:func:`matEleSelector` function.


    Examples
    --------
    # Calculate for standard matE dataset, restacked for [l,m,mu].
    >>> daOut, daDot = densityCalc(matE, denDims = ['LM', 'mu'], selDims = {'Type':'L'}, thres = 1e-1)
    # Calculate for standard matE, no restack and sum all other dims.
    >>> daOut, daDot = densityCalc(matE, denDims = 'LM', selDims = {'Type':'L'}, thres = 1e-1, sumDims = True)
    # Calculate for standard matE, no restack and sum all other dims, except Sym.
    >>> daOut, daDot = densityCalc(matE, denDims = 'LM', selDims = {'Type':'L'}, thres = 1e-1, sumDims = True, keepDims = 'Sym')

    """

    # Set data
    denSettings = locals()
    attrs = da.attrs.copy()
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

    # Set rsMap for singleton dim cases (not set in dimRestack, but maybe should be)
    if not rsMap:
        if dimCheck['shared']:
            rsMap = {denDim:[denDim]}  # Case for single, unstacked dim set
        else:
            rsMap = {denDim:dimCheck['stackedMap'][denDim]}  # For case of single, already stacked dim, dimCheck returns rsMap = {}.

    newDims = {item:item+'_p' for item in rsMap[denDim]}
    denDimP = denDim+'_p'
    denDimPMap = {denDimP:list(newDims.values())}
    daConj = daDen.conj().unstack(denDim).rename(newDims).stack(denDimPMap)
    daDot = daDen * daConj

    #*** Propagate attrs & set outputs
    daDot.attrs = attrs

    daDot.attrs['dataType']='Density Matrix'
    daDot.attrs['density'] = {'denSettings':denSettings,
                             'sumTot':sumTot,
                             'denDims':rsMap.update(denDimPMap),
                             'kdims':[denDim, denDimP]}

    # Note denSettings will include input da, can skip with, e.g.:
#     [BetasNormX.attrs.update({k:calcSettings[k]}) for k in calcSettings.keys() if not (isinstance(calcSettings[k], xr.DataArray))]  # Slightly ugly, but set only items which are not Xarrays.

    daDot.attrs['kdims']=daDot.attrs['density']['kdims']

    # General .dot version... might be faster than above?
#     matEdot = xr.dot(daDen, daDen.conj().rename({denDim:denDim+'_p'}), dims = sumDims)

#     if sumTot:
#         daDot = daDot.sum(sumTot, keep_attrs = True,  skipna = ~keepNaNs)  # Only run if sumTot not empty?
                                    # Otherwise .sum([]) will push NaNs > zeros

    daOut = matEleSelector(daDot.sum(sumTot, keep_attrs = True,  skipna = ~keepNaNs),
                           thres = thres, sq = squeeze) # Threshold density mat
                                                                        # TODO: pass **kwargs here?
                                                                        # Pass dims = denDims?

#     daOut.attrs = daDot.attrs  # Propagate attrs - may be lost after .sum(), should be OK if .sum(keep_attrs = True)

    return daOut, daDot


#******* HV plotting code
#
#  TODO:
#  - Move to hvPlotters
#  - Add heatmap defaults to hvPlotters
#  - Seaborn version for matrix plotting.

import xarray as xr
import pandas as pd
import holoviews as hv
# import hvplot.pandas
hv.extension('bokeh')

from epsproc import plotTypeSelector

# Quick default settings from tmo-dev, tmoDataBase.py init.

imgSize = 600
from holoviews import opts

opts.defaults(opts.HeatMap(width=imgSize, frame_width=imgSize, aspect='square', tools=['hover'], colorbar=True))



def matPlot(da, kdims = None, pTypes = ['r','i','a'],
            negQs = True, thres = None,
            returnType = 'plot', verbose = 1):
    """
    General matrix (2D) plot + stacked dims plotter with HoloViews.

    Parameters
    ----------

    da : Xarray
        Input dataset for plotting.
        Assumed to contain a "matrix" (2D) to plot, with other dims stacked to HV plot widgets.

    kdims : list, optional, default = None
        Key dims for plotting.
        If None,
        - will be taken from da.attrs[kdims] if it is set.
        - otherwise, uses da.dims[-2:]

    pTypes : list, optional, default = ['r','i','a']
        List of plot types to set (via :py:func:`epsproc.plotTypeSelector`_
        These will be available via the HV plot widgets.

    negQs : NOT YET IMPLEMENTED
        Include -ve valued QNs in plot?

    thres : NOT YET IMPLEMENTED
        optional, float, default = None
        Threshold for plot.

    returnType : str, optional, default = 'plot'
        If 'plot' return plot object only.
        Otherwise, returns plot + datasets = (hvmap, hvds, daPlotDS)

    verbose : int, bool, optional, default = 1
        Print additional info if true.

    Returns
    -------
    Holoviews holomap object
        If 'plot'; otherwise, returns plot + datasets = (hvmap, hvds, daPlotDS)


    """

    #*** Set data
    daPlot = da.copy()  # May want to add thresholds etc. here?
    attrs = daPlot.attrs.copy()

    if not daPlot.name:
        daPlot.name = 'Unnamed'  # Set default name if blank, otherwise HV will error.

    # Try and set kdims if not passed
    if kdims is None:
        if 'kdims' in attrs.keys():
            kdims = attrs['kdims']
        else:
            kdims = list(daPlot.dims[-2:])  # Use final two dims? This matches denCalc outputs.

        if verbose:
            print(f'Set plot kdims to {kdims}; pass kdims = [dim1,dim2] for more control.')


#     if not negQs:
#         if daPlot.dims


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
    daPlotDS = xr.Dataset({pType:plotTypeSelector(daPlot.copy(), pType = pType) for pType in pTypes}) #['r','i','a']})
    daPlot = daPlotDS.to_array().rename({'variable':'pType'})  # Restack pType to array
    daPlot.attrs = attrs  # Propagate attrs
    daPlot.name = da.name
    hvds = hv.Dataset(daPlot)
    hvmap = hvds.to(hv.HeatMap, kdims=kdims)

    if returnType is 'full':
        return hvmap, hvds, daPlotDS
    else:
        return hvmap
