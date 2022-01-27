"""
ePSproc basic plotting functions

Some basic functions for 2D/3D plots.

19/11/19        Adding plotting routines for matrix elements & beta values, to supercede existing basic methods.

07/11/19    v1  Molecular plotter for job info.


TODO
----


"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# For Arrow3D class
from matplotlib.patches import FancyArrowPatch
from mpl_toolkits.mplot3d import proj3d

# Package functions
from epsproc.util.listFuncs import BLMdimList
from epsproc.sphPlot import plotTypeSelector
# from epsproc.util import matEleSelector  # Throws error, due to circular import refs?
# from epsproc.util.conversion import multiDimXrToPD
# from epsproc.util import matEdimList, BLMdimList, dataTypesList, multiDimXrToPD

# Additional plotters
try:
    import seaborn as sns
    import epsproc._sns_matrixMod as snsMatMod  # SNS code with modified clustermap - NOTE this may need Seaborn 0.9.0 - TBC
except ImportError as e:
    if e.msg != "No module named 'seaborn'":
        raise
    print('* Seaborn not found, clustermap plots not available. ')

# Holoviews backend routines - already has appropriate try/except in code
# TODO: may want to move more setup routines to .plot?
import epsproc.plot.hvPlotters as hvPlotters
from epsproc.plot.util import showPlot
hvFlag = hvPlotters.hvFlag
setPlotters = hvPlotters.setPlotters

print('* Setting plotter defaults with epsproc.basicPlotters.setPlotters(). Run directly to modify, or change options in local env.')
setPlotters()

# if hvFlag:
#     setPlotters()
#     print()

# hv = hvPlotters.hv


# Env check
from epsproc.util.env import isnotebook
__notebook__ = isnotebook()

# Arrow3D class from https://stackoverflow.com/questions/22867620/putting-arrowheads-on-vectors-in-matplotlibs-3d-plot
# Code: https://stackoverflow.com/a/22867877
# Builds on existing 2D arrow class, projects into 3D.
# From StackOverflow user CT Zhu, https://stackoverflow.com/users/2487184/ct-zhu
class Arrow3D(FancyArrowPatch):
    """
    Define Arrow3D plotting class

    Code verbatim from StackOverflow post https://stackoverflow.com/a/22867877
    Thanks to CT Zhu, https://stackoverflow.com/users/2487184/ct-zhu

    """

    def __init__(self, xs, ys, zs, *args, **kwargs):
        FancyArrowPatch.__init__(self, (0,0), (0,0), *args, **kwargs)
        self._verts3d = xs, ys, zs

    def draw(self, renderer):
        xs3d, ys3d, zs3d = self._verts3d
        xs, ys, zs = proj3d.proj_transform(xs3d, ys3d, zs3d, renderer.M)
        self.set_positions((xs[0],ys[0]),(xs[1],ys[1]))
        FancyArrowPatch.draw(self, renderer)


# Basic plotting for molecular structure
# TODO: replace with more sophisticated methods!
def molPlot(molInfo):
    """Basic 3D scatter plot from molInfo data."""

    #*** Plot atoms
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    scatter = ax.scatter(molInfo['atomTable'][:,2], molInfo['atomTable'][:,3], molInfo['atomTable'][:,4], s = 100* molInfo['atomTable'][:,0], c = molInfo['atomTable'][:,0])
    ax.axis('off');

    # Add labels + bonds for small systems (need more logic for larger systems!)
    if molInfo['atomTable'].shape[0] < 4:
        ax.plot(molInfo['atomTable'][:,2], molInfo['atomTable'][:,3], molInfo['atomTable'][:,4])

        lShift = 0.01* np.array([1., 1., 1.])
        for row in molInfo['atomTable']:
            # ax.text(row[2] + lShift[0], row[3] + lShift[1], row[4] + lShift[2], f"Z = {row[0]}") # molInfo['atomTable'][:,0])
            ax.text(row[2] + lShift[0], row[3] + lShift[1], row[4] + lShift[2], f"Z = {row[0].data}") # Version for Xarray input

    else:
        # Generate legend from scatter data, as per https://matplotlib.org/3.1.1/gallery/lines_bars_and_markers/scatter_with_legend.html
        legend1 = ax.legend(*scatter.legend_elements(), loc="lower left", title="Z")
        ax.add_artist(legend1)


    #*** Plot axes with direction vecs

    dirs = ('x','y','z')  # Labels
    oV = np.zeros([3,3])   # Generate origin
    xyzV = np.identity(3)  # Generate unit vectors

    # Scale for lines, check axis limits & molecular extent
    # sf = 1
    axLims = np.asarray([ax.get_xlim(), ax.get_ylim(), ax.get_zlim()])
    sfAxis = np.asarray([ax.get_xlim()[1], ax.get_ylim()[1], ax.get_zlim()[1]])  # Get +ve axis limits
    sfMol = 1.8* np.asarray([molInfo['atomTable'][:,2].max(), molInfo['atomTable'][:,3].max(), molInfo['atomTable'][:,4].max()]) # Get max atom position

    # Set sensible limits
    mask = sfMol > 0
    sfPlot = sfMol*mask + (sfAxis/2)*(~mask)

    for n, xyz in enumerate(xyzV):
        sf = sfPlot[n]
    #     ax.plot([oV[n,0], sf*xyz[0]], [oV[n,1], sf*xyz[1]], [oV[n,2], sf*xyz[2]])  # Basic line
    #     ax.quiver(oV[n,0],oV[n,1],oV[n,2],sf*xyz[0],sf*xyz[1],sf*xyz[2], arrow_length_ratio = 0.05, lw=2)  # Use quiver for arrows - not very pleasing as arrow heads suck
        a = Arrow3D([oV[n,0], sf*xyz[0]], [oV[n,1], sf*xyz[1]], [oV[n,2], sf*xyz[2]], mutation_scale=10,
                    lw=4, arrowstyle="-|>", color="r", alpha=0.5)  # With Arrow3D, defined above
        ax.add_artist(a)
        ax.text(sf*xyz[0], sf*xyz[1], sf*xyz[2], dirs[n])


    # Add (y,z) plane as a surface plot
    yy, zz = np.meshgrid(range(-1, 2), range(-1, 2))
    xx = np.zeros([3,3])
    ax.plot_surface(xx*sfPlot[0],yy*sfPlot[1],zz*sfPlot[2], alpha=0.1)

    # Reset plot limits (otherwise rescales to arrows!)
    ax.set_xlim(axLims[0,:])
    ax.set_ylim(axLims[1,:])
    ax.set_zlim(axLims[2,:])

    # TO CONSIDER:
    # See https://matplotlib.org/mpl_toolkits/mplot3d/api.html#mpl_toolkits.mplot3d.axes3d.Axes3D
    # ax.get_proj()  # Get projection matrix
    # ax.view_init(0,20)  # Set view, (az,el)

    plt.show()



#*** BLM surface plots

def BLMplot(BLM, thres = 1e-2, thresType = 'abs', selDims = None,
            XS = False, xDim = 'Eke', pType = 'r',
            col = 'Labels', row = None, pStyle = 'line',
            backend = 'xr', returnPlot = False, renderPlot = True,
            **kwargs):
    """
    Plotting routines for BLM values from Xarray.
    Plot line or surface plot, with various backends available.

    Parameters
    ----------
    BLM : Xarray
        Input data for plotting, dims assumed to be as per :py:func:`ePSproc.util.BLMdimList()`

    thres : float, optional, default 1e-2
        Value used for thresholding results, only values > thres will be included in the plot.
        Either abs or relative (%) value, according to thresType setting.
        Set thres = None to skip.

    thresType : str, optional, default = 'abs'
        Set to 'abs' or 'pc' for absolute threshold, or relative value (%age of max value in dataset)

    selDims : dict, optional, default = None
        Dimensions to select from input Xarray.

    XS : bool, optional, default = False
        Use absolute XS value for l=0 if true (from BLM.XS); otherwise show current values (usually normalised to =1).

    xDim : str, optional, default = 'Eke'
        Dimension to use for x-axis, also used for thresholding. Default plots (Eke, BLM) surfaces with subplots for (Euler).
        Change to 'Euler' to plot (Euler, BLM) with (Eke) subplots.

    col, row : strings, optional, default = 'Labels', None
        Set for dim handling in plot.
        For xr backend these define plot layout.

    pType : char, optional, default 'a' (abs values)
        Set (data) type to plot. See :py:func:`plotTypeSelector`.

    pStyle : string, optional, default = 'line'
        NOT YET PROPERLY IMPLEMENTED
        For XR backend can also set 'surf'.

    backend : str, optional, default = 'xr'
        Plotter to use. Currently supports 'xr' = Xarray, 'hv' = Holoviews
        Default is 'xr' for Xarray internal plotting. May be switched according to plot type in future...
        **kwargs are currently passed to the plotter backend.

    returnXR : bool, optional, default = False
        Return plot data to caller if true.

    returnPlot : bool, optional, default = False
        For hv backend, return object to caller?

    renderPlot : bool, optional, default = True
        For hv backend, render plot (if running in a notebook)?

    Notes
    -----
    Updates in issue #27: https://github.com/phockett/ePSproc/issues/27
    Demo: https://epsproc.readthedocs.io/en/dev/tests/hvPlotters_fn_tests_150720_v250122.html

    - Proper dim handling to be implemented. See code elsewhere(?).
        XR plotter currently uses 'cDims' in some places too.
    - Add XR and HV data return, currently only return objects.

    27/01/22: Updated with better selection & basic XR dim handling, but still needs work.
              Added basic HV plotter routines.

    """
    # Local/deferred import to avoid circular import issues at module level.
    # TO DO: Should fix with better __init__!
    from epsproc.util import matEleSelector

    # Check dims, and set facet dim
    # TODO: update with newer dim handling routines (dynamic)
    dims = BLMdimList()
    dims.remove(xDim)
    cDims = dims.remove('BLM')

    # For %age case
    if thresType is 'pc':
        thres = thres * BLM.max()

    # Threshold & set data for plotting. Set to check along Eke dim, but may want to pass this as an option.
    BLMplot = matEleSelector(BLM, thres=thres, inds = selDims, dims = xDim, sq = True)

    # Set dataType if missing
    if 'dataType' not in BLMplot.attrs:
        BLMplot.attrs['dataType'] = '(No dataType)'
        # print(f"Set dataType {daPlot.attrs['dataType']}")

    if XS and (BLMplot.attrs['dataType'] is 'BLM'):
        # daPlot.values = daPlot * daPlot.XS
        BLMplot = BLMplot.where(BLMplot.l !=0, BLMplot.XS)  # Replace l=0 values with XS

    BLMplot = plotTypeSelector(BLMplot, pType = pType, axisUW = xDim)

    #*** Fix any missing attrs - set this LAST otherwise may be accidentally overwritten above.
    # Set filename if missing
    if 'file' not in BLMplot.attrs:
        BLMplot.attrs['file'] = '(No filename)'

    # For HV plotters need a name!
    # if hasattr(BLMplot,'name') and not BLMplot.name:
    if not BLMplot.name:
        BLMplot.name = 'BLM'

    #*** Plot
    if backend is 'xr':

        # Set BLM index for plotting. TODO: check for updated XR versions.
        # Necessary for Xarray plotter with multi-index categories it seems.
        BLMplot['BLMind'] = ('BLM', np.arange(0, BLMplot.BLM.size))

        if pStyle is 'surf':
            # BLMplot.real.squeeze().plot(x=xDim, y='BLMind', col=cDims, size = 5, **kwargs)
            BLMplot.plot.line(x=xDim, y='BLMind', col=cDims, row=row, **kwargs)
            # .real.plot.line(x=Etype, col=col, row=row, **kwargs)

        else:
            # BLMplot.plot.line(x=xDim, col=cDims, row=row, **kwargs)
            BLMplot.plot.line(x=xDim, col=col, row=row, **kwargs)

    if backend is 'hv':
        # Basic HVplotter routine, need to add better dim handling!
        # hvObj = hvPlotters.curvePlot(BLMplot.unstack('BLM').squeeze().fillna(0), kdims='Eke', renderPlot = False)
        hvObj = hvPlotters.curvePlot(BLMplot.unstack('BLM'), kdims='Eke', renderPlot = False)  # OK, don't render

        # Use showPlot() to control render & return, or just return object - might be overkill?
        # Note default .overlay(['l','m'])
        if renderPlot:
            return showPlot(hvObj.overlay(['l','m']), returnPlot = returnPlot, __notebook__ = __notebook__)  # Currently need to pass __notebook__?

        else:
            return hvObj


#************************* Routines for matrix element plotting

# Function for getting list of unique syms
def symListGen(data):
    symList = []
    [symList.append(list(item)) for item in data.get_index('Sym').to_list()]

    return np.unique(np.ravel(symList))

def lmPlot(data, pType = 'a', thres = 1e-2, thresType = 'abs', SFflag = True, logFlag = False, eulerGroup = True,
        selDims = None, sumDims = None, plotDims = None, squeeze = False, fillna = False,
        xDim = 'Eke', backend = 'sns', cmap = None, figsize = None, verbose = False, mMax = 10, titleString = None,
        labelRound = 3):
    """
    Plotting routine for ePS matrix elements & BLMs.

    First pass - based on new codes + util functions from sphPlot.py, and matE sorting codes.

    Parameters
    -----------
    data : Xarray, data to plot.
        Should work for any Xarray, but optimised for dataTypes:
        - matE, matrix elements
        - BLM paramters
        - ADMs

    pType : char, optional, default 'a' (abs values)
        Set (data) type to plot. See :py:func:`plotTypeSelector`.

    thres : float, optional, default 1e-2
        Value used for thresholding results, only values > thres will be included in the plot.
        Either abs or relative (%) value, according to thresType setting.

    thresType : str, optional, default = 'abs'
        Set to 'abs' or 'pc' for absolute threshold, or relative value (%age of max value in dataset)

    SFflag : bool, optional, default = True
        For dataType = matE: Multiply by E-dependent scale factor.
        For dataType = BLM: Multiply by cross-section (XS) (I.e. if False, normalised BLMs are plotted, with B00 = 1)

    logFlag : bool, optional, default = False
        Plot values on log10 scale.

    eulerGroup : bool, optional, default = True
        Group Euler angles by set and use labels (currenly a bit flakey...)

    selDims : dict, optional, default = {'Type':'L'}
        Dimensions to select from input Xarray.

    sumDims : tuple, optional, default = None
        Dimensions to sum over from the input Xarray.

    plotDims : tuple, optional, default = None
        Dimensions to stack for plotting, also controls order of stacking (hence sorting and plotting).
        Default case plots all dims not included in (xDim, sumDims, selDims).
        E.g. for matrix elements this  will be ('l','m','mu','Cont','Targ','Total','it','Type')
        TO DO: auto generation for different dataType, also based on selDims and sumDims selections.
        11/03/20: partially fixed, now add any missing dims to plot automatically.
        NOTE: this currently fails for 'Labels' case, not sure why.

    squeeze : bool, optional, default = True
        Drop singleton dimensions from plot.

    fillna : bool, optional, default = False
        Fill NaN values with 0 if True.

    xDim : str, optional, default = 'Eke'
        Dimension to use for x-axis, also used for thresholding. Default plots (Eke, LM) surfaces.
        NOTE: if xDim is a MultiIndex, pass as a dictionary mapping, otherwise it may be unstacked during data prep.
        E.g. for plotting stacked (L,M), set xDim = {'LM':['L','M']}
        NOTE: this is currently also passed to matEleSelector(), so stacked dims *must* exist in inital Xarray.
        THIS SHOULD BE FIXED, since it's a pain.

    backend : str, optional, default = 'sns'
        Plotter to use. Default is 'sns' for Seaborn clustermap plot.
        Set to 'xr' for Xarray internal plotting (not all passed args will be used in this case). May be switched according to plot type in future...

    cmap : str, optional, default = None
        Cmap option to pass to sns clustermap plot.

    figsize : tuple, optional, default None
        Tuple for Seaborn figure size (ratio), e.g. figsize = (15,5).
        Useful for setting a long axis explicitly in cases with large dimensional disparity.
        Default results in a square (ish) aspect.

    verbose : bool, optional, default False
        Print debug info.

    mMax : int, optional, default = 10
        Used as default for m colour mapping, in cases where coords are not parsed.
        TODO: fix this, it's ugly.

    titleString : str, optional, default = None
        Set main title for plot. If not set, filename or data.attr['jobLabel'] will be used if present.
        Some plot settings labels will also be appended to this.

    labelRound : int, optional, default = 3
        Round to N d.p. for plot labels (legend text).
        This is also passed to :py:func:`epsproc.multiDimXrToPD()` for column labels (for E, t or Euler dims only)

    Returns
    -------
    daPlot : Xarray
        Data subset as plotted.

    legendList : list
        Labels & colour maps

    g : figure object


    Notes
    -----
    * Data is automagically sorted by dims in order set in plotDim.
    * For clustermap use local version - code from Seaborn, version from PR1393 with Cluster plot fixes.
        * https://github.com/mwaskom/seaborn/pull/1393
        * https://github.com/mwaskom/seaborn/blob/fb1f87e800e69ba2e9309f922f9dac470e3a6c78/seaborn/matrix.py
        * This may break for versions of Seaborn other than 0.9.0 (tested version)
    * Currently only set for single colourmap choice, should set as dict.
    * Clustermap methods from:
        * https://stackoverflow.com/questions/27988846/how-to-express-classes-on-the-axis-of-a-heatmap-in-seaborn
        * https://seaborn.pydata.org/examples/structured_heatmap.html
    * Seaborn global settings currently included here:
        * sns.set(rc={'figure.dpi':(120)})
        * sns.set_context("paper")
        * These are reset at end of routine, apart from dpi.

    To do
    -----
    - Improved dim handling, maybe use :py:func:`epsproc.util.matEdimList()` (and related functions) to avoid hard-coding multiple cases here.
        - Partially implemented in dimMap.
        - Currently throws an error for None in symmetry plotting cases for unrecognised dataType. TO FIX!
           - 14/01/21 fixed, but possibly introduced other issues!
        - Consider consolidating mapping styles further, and keeping a global list of parent and child dims? Or otherwise pulling from stacked Xarray?
    - Improved handling of sets of polarization geometries (angles).
    - CONSOLIDATE stacked/unstacked dim handling.  At the moment some functions use stacked, some unstacked, which is leading to annoying bugs.

    History
    -------

    * 14/01/21 - Fixed bug in handling nested dims for unknown datatypes (hopefully), specifically for Sym dims.
    * 27/02/20 - handling for MultiIndex cases for Wigner 3j type data (function of (l,lp,L,m,mp,M))
        * Fixed issue with conversion to Pandas table - should now handle dims better.
        * Improved multidim dropna() using Xarray version (previously used Pandas version).
        * Added xDim as dict mapping capability.
    * 05/12/19 - debug
    * 29/11/19 - initial version.

    Examples
    --------
    (See https://github.com/phockett/ePSproc/blob/master/epsproc/tests/ePSproc_demo_matE_plotting_Nov2019.ipynb)

    """
    # Local/deferred import to avoid circular import issues at module level.
    # TO DO: Should fix with better __init__!
    from epsproc.util import matEleSelector, matEdimList, BLMdimList, dataTypesList, multiDimXrToPD

    # Set Seaborn style
    # TO DO: should pass args here or set globally.
    # sns.set()
    sns.set_context("paper")  # "paper", "talk", "poster", sets relative scale of elements
                            # https://seaborn.pydata.org/tutorial/aesthetics.html
    # sns.set(rc={'figure.figsize':(11.7,8.27)})  # Set figure size explicitly (inch)
                            # https://stackoverflow.com/questions/31594549/how-do-i-change-the-figure-size-for-a-seaborn-plot
                            # Wraps Matplotlib rcParams, https://matplotlib.org/tutorials/introductory/customizing.html
    sns.set(rc={'figure.dpi':(120)})

#*** Data prep
    # Make explicit copy of data
    if selDims is not None:
        daPlot = data.sel(selDims).copy()   # Select here instead of via matEleSelector, to avoid issues with %age threshold case.
    else:
        daPlot = data.copy()

    daPlot.attrs = data.attrs
#     daPlot = data

    # Set dataType if missing - NOTE additionally set at end to avoid accidentally dropping this.
    if 'dataType' not in daPlot.attrs:
        daPlot.attrs['dataType'] = '(No dataType)'
        print(f"Set dataType {daPlot.attrs['dataType']}")

    # Use SF (scale factor)
    # Write to data.values to make sure attribs are maintained.
    if SFflag and (daPlot.attrs['dataType'] is 'matE'):
        daPlot.values = daPlot * daPlot.SF

    if SFflag and (daPlot.attrs['dataType'] is 'BLM'):
        daPlot.values = daPlot * daPlot.XS

    # For %age case
    if thresType is 'pc':
        thres = thres * np.abs(daPlot.max()).values  # Take abs here to ensure thres remains real (float)
                                              # However, does work without this for xr.where() - supports complex comparison.

    # Get full dim list
    # Manual
    # if daPlot.attrs['dataType'] is 'matE':
    #     dimMap = matEdimList(sType = 'sDict')
    # else:
    #     dimMap = BLMdimList(sType = 'sDict')

    # Using util.dataTypesList() - THIS IS NOT yet used in the main plotting routine.
    # ACTUALLY - is used for Sym plotting case, may throw an error for None case.
    # See lines ~630, 696
    dataTypes = dataTypesList()
    if daPlot.attrs['dataType'] in dataTypes:
        dimMap = dataTypes[daPlot.attrs['dataType']]['dims']
    else:
        # dimMap = None  # This is painful and gives errors for arb data types symmetry cases
        # dimMap = list(daPlot.dims)  # This is OK, but won't get nested dims.
        dimMap = {item:(list(daPlot.get_index(item).names) if len(daPlot.get_index(item).names)>1 else []) for item in daPlot.dims }  # Return subdims or empty list in dict format



    # Eulers >>> Labels
    if eulerGroup and ('Euler' in daPlot.dims):
        # if 'Labels' in daPlot.drop('Euler').dims:  # This did work, but no longer (Jan 2020) - may be an Xarray version issue?  Using cords should be safe however.
        if 'Labels' in daPlot.coords:
            daPlot = daPlot.drop('Euler').swap_dims({'Euler':'Labels'})   # Set Euler dim to labels
        else:
            pass
            # TODO: add labels here if missing. See setPolGeoms()
            # Set labels if missing, alphabetic or numeric
            # if eulerAngs.shape[0] < 27:
            #     labels = list(string.ascii_uppercase[0:eulerAngs.shape[0]])
            # else:
            #     labels = np.arange(1,eulerAngs.shape[0]+1)



# Restack code from mfblm()
    # # Unstack & sub-select data array
    # daUnStack = da.unstack()
    # daUnStack = matEleSelector(daUnStack, thres = thres, inds = selDims, sq = True)
    #
    # # Check dims vs. inds selection and/or key dims in output
    # for indTest in sumDims:
    #     # if (indTest in inds) or (indTest not in unstackTest.dims):
    #     if indTest not in daUnStack.dims:
    #         daUnStack = daUnStack.expand_dims(indTest)
    #
    # # Restack array along sumInds dimensions, drop NaNs and squeeze.
    # daSumDim = daUnStack.stack(SumDim = sumDims).dropna('SumDim').squeeze()

    # Sum & threshold
    if sumDims is not None:
        daPlot = daPlot.sum(sumDims) #.squeeze()  # Set squeeze condition below.
        daPlot.attrs = data.attrs  # Reset attribs

    # If xDim is stacked, check it exists and create it if not.
    # TOTAL FUDGE HERE - this may be problematic due to unstack() UGH.
    # HOWEVER, without this xDim must exist in passed daPlot array, otherwise later code may fail.
    if type(xDim) == dict:
        dimCheck = any([item in xDim for item in daPlot.dims])
        if not dimCheck:
            daPlot = daPlot.unstack().stack(xDim)

    # Threshold on abs() value before setting type, otherwise all terms will appear for some cases (e.g. phase plot)
    # daPlot = matEleSelector(daPlot, thres=thres, inds = selDims, dims = xDim) # , sq = True)  # Squeeze may cause issues here if a singleton dim is used for xDim.
    daPlot = matEleSelector(daPlot, thres=thres, dims = xDim, drop = False, sq = squeeze) # Version without selection, now incorporated above.
    daPlot = plotTypeSelector(daPlot, pType = pType, axisUW = xDim)

    # daPlot = ep.util.matEleSelector(daPlot, thres=thres, inds = selDims, dims = 'Eke', sq = True)

    if logFlag:
        daPlot.values = daPlot.pipe(np.log10)

#*** Fix any missing attrs - set this LAST otherwise may be accidentally overwritten above.
    # Set filename if missing
    if 'file' not in daPlot.attrs:
        daPlot.attrs['file'] = '(No filename)'

    # Set dataType if missing
    if 'dataType' not in daPlot.attrs:
        daPlot.attrs['dataType'] = '(No dataType)'
        print(f"Set dataType {daPlot.attrs['dataType']}")


#*** Plotting routines
    # Rough code for xr plotter.
    if backend is 'xr':
        # Set index for plotting
        daPlot['LMind'] = ('LM',np.arange(0, daPlot.LM.size))

        # Plot abs values, with faceting on symmetry (all mu)
        # TODO: faceting with passed or smart dims.
        daPlot.plot(x='Eke', y='LMind', col='Sym', row='Type', robust = True)


    if backend is 'sns':

        print(f"Plotting data {daPlot.attrs['file']}, pType={pType}, thres={thres}, with Seaborn")

        # *** Dims & cmaps
        # Get sym list if applicable, this is used for unique colour-mapping over dims and legends.
        if 'Sym' in daPlot.dims:
            symList = symListGen(daPlot)
            symFlag = True
            palSym = sns.color_palette("hls", np.unique(symList).size)  # Generate unified cmap
            lutSym = dict(zip(map(str, np.unique(symList)), palSym))
        else:
            symFlag = False
            symList = []

        # 14/01/20 Set catch-all defaults for dims with shared cmaps - UGLY, partially implemented below, but still have lots of special cases
        # 01/02/21 Added set check with xDim to ensure this works for mixed (l,K) cases, with one as xDim (otherwise can miss labels).
        lDims = list(set(['l', 'K', 'lp', 'L'])-set(xDim))
        mDims = ['m', 'mp', 'mu', 'Q', 'S', 'mu0', 'Lambda', 'M']

        # Set unified cmap for (m,mu)
        # Switch for matE or BLM data type.
#         if ('m' in daPlot.dims) and ('mu' in daPlot.dims):

        # Set a default case to allow for arb. dim names
        # TODO: fix this to select based on max in passed QNs.
        mList = np.arange(-mMax, mMax+1)

        # Set unstaked (full) dim list
        dimUS = daPlot.unstack().dims

        # Set l (or equivalent) cmap here, ensures ordering over over dims (ordered list returned from daPlot Xarray).
        # NOTE - set crude default case for max = 10.
        # NOTE - assumes only one of these dims is present (or will overwrite)
        # Added 28/09/20
        # lList = []
        lList = np.arange(0, 10, dtype=np.int64)  # Set default case. UGLY.
        # lList = daPlot[xDim].pipe(np.unique)        # Should be able to default to first xDim passed?
        lListTemp = []
        for dim in dimUS:
            # if dim in ['l', 'K', 'lp', 'L']:
            if dim in lDims:
                # lList = daPlot[dim].pipe(np.unique)  # This fails for multiple l-types with different values
                lListTemp.extend(daPlot[dim].pipe(np.unique))  # Stack over all matching dims then recheck for unique values.

        if lListTemp:
            lList = np.unique(lListTemp)

        # if not lList.any():
        #     lList = lListDefault


        palL = sns.light_palette("green", n_colors=lList.size)
        lutL = dict(zip(map(str, lList), palL))


        if 'LM' in daPlot.dims:
            try:
                mList = np.unique(daPlot.LM.m)   # For Wigner 3j functions, allow for m or M
            except AttributeError as e:
                pass
                # if e != "'DataArray' object has no attribute 'm'":
                #     raise
            try:
                mList = np.unique(daPlot.LM.M)
            except AttributeError as e:
                # if e != "'DataArray' object has no attribute 'M'":
                #     raise
                pass

        if 'BLM' in daPlot.dims:
            mList = np.unique(daPlot.BLM.m)

        if 'ADM' in daPlot.dims:
            mList = np.unique([daPlot.ADM.Q, daPlot.ADM.S])

        mColours = mList.size

        if mColours < 3:  # Minimal value to ensure mu mapped - may be better to generate without ref here?
            mColours = 3

#         palM = sns.diverging_palette(220, 20, n=mColours) # sns.color_palette("hls", mList.size)
        palM = sns.color_palette("coolwarm", mColours)
        lutM = dict(zip(map(str, mList), palM))

# 12/03/20 - moving this to separate conversion function, since tabulation is handy.
        # Note thres=None set here to avoid putting holes in the pd data. This is a temp workaround!
        daPlotpd, daPlot = multiDimXrToPD(daPlot, colDims = xDim, rowDims = plotDims, thres = thres, dropna = True, fillna = fillna, squeeze = squeeze, colRound = labelRound)

        # # Convert to Pandas 2D array, stacked along plotDims - these will be used for colour bar.
        # # TO DO: fix hard-coded axis number for dropna()
        # # 27/02/20: switched to use xr.dropna(dim = 'plotDim', how = 'all') instead of pd.dropna.  This should be more robust... hopefully won't break old code.
        # #           Also added stack(xDim = xDim) to allow for multilevel X plotting.
        # # 11/03/20: added auto setting for plotDims, based on sets, plus checks for missing dims.
        #
        # # Work-around for dict case - get list of unstacked dims for comparison
        # if type(xDim) == dict:
        #     xDimList = list(xDim.items())[0][1]
        # else:
        #     xDimList = xDim
        #
        # if plotDims is None:
        #     plotDims = list(set(dimUS) - set(xDimList))   # Use set arithmetic to get items
        #     plotDims.sort()                      # Set sort to return alphebetical list.
        #
        #
        # # Check plotDims exist, otherwise may throw errors with defaults
        # plotDimsRed = []
        # # for dim in daPlot.unstack().dims:
        # #     if dim in plotDims:
        # #         plotDimsRed.append(dim)
        # for dim in plotDims:
        #     if dim in dimUS:
        #         plotDimsRed.append(dim)
        #
        # # Additional check for any missing dims
        # # Check # of dims and correct for any additional/skipped dims
        # # Bit ugly - should be integrated with above code
        # if (len(xDimList) + len(plotDimsRed)) != len(dimUS):
        #     for dim in dimUS:
        #         if not (dim in xDimList) and not (dim in plotDimsRed):
        #             plotDimsRed.append(dim)
        #
        #             if verbose:
        #                 print(f'Adding {dim} to plotting dim list.')
        #
        #
        # # Restack for plotting, and drop singleton dimensions if desired.
        # daPlot = daPlot.unstack().stack(plotDim = plotDimsRed).dropna(dim = 'plotDim', how = 'all')
        #
        # # Restack xDim in cases where it is a MultiIndex
        # if type(xDim) == dict:
        #     daPlot = daPlot.stack(xDim)
        #
        #
        # if squeeze:
        #     # daPlotpd = daPlot.unstack().stack(plotDim = plotDimsRed).squeeze().to_pandas().dropna(axis = 1).T
        #     # daPlotpd = daPlot.unstack().stack(plotDim = plotDimsRed).dropna(dim = 'plotDim', how = 'all').squeeze().to_pandas().T
        #     daPlotpd = daPlot.squeeze().to_pandas()
        #
        # else:
        #     # daPlotpd = daPlot.unstack().stack(plotDim = plotDimsRed).to_pandas().dropna(axis = 1).T
        #     # daPlotpd = daPlot.unstack().stack(plotDim = plotDimsRed).dropna(dim = 'plotDim', how = 'all').to_pandas().T
        #     daPlotpd = daPlot.to_pandas()
        #
        # # Transpose Pandas table if necessary - xDim must be columns
        # if type(xDim) != dict:
        #     if xDim not in daPlotpd.columns.names:
        #         daPlotpd = daPlotpd.T


        # For dictionary case, check items for each key are in column names.
        # THIS CODE IS HORRIBLE - should be a neater way to do this.
        # TODO: fix case for some levels missing, at the moment assumes any present is OK.
        # TODO: test for single vs. MultiIndex case - columns.names vs. columns.name?
        # else:
        #     for key in xDim:
        #         dimList = xDim[key]
        #         check = [item in daPlotpd.columns.names for item in dimList]
        #         if not any(check):
        #             daPlotpd = daPlotpd.T
# --- END pandas conversion.

        # For Wigner3j datatypes, still get some illegal values creeping through - now try removing again...
        # if daPlot.dataType == 'Wigner3j':
        #     daPlotpd = daPlotpd.dropna(axis = 1, how = 'all')
        # NOW set in multiDimXrToPD, with dropna = True (default)

        # Set multi-index indicies & colour mapping
        cList = []
        legendList = []
        for dim in daPlotpd.index.names:
            # 28/09/20: added .unique().sort_values(), should clean up legend a little bit?
            # EDIT: actually added only for legendList.append() stuff below.
            # Q: this should be OK, but may mess up multidim labels...?
            Labels = daPlotpd.index.get_level_values(dim)  # .unique().sort_values() #.astype('str')

            if verbose:
                print(dim)
                print(Labels.unique().size)

            # For (l,m,mu) set colour scales as linear (l), and diverging (m,mu)
            # May want to modify (m,mu) mappings to be shared?
            # 28/09/20: revisiting to fix ordering over multiple states.
            #   Should be fixed by setting master list above, as per M state settings, since that takes ordered list over full dataset.
            #   (Otherwise l cmap may be split by sym or other dim ordering)
            # if dim in ['l', 'K', 'lp', 'L']:
            if dim in lDims:
                # pal = sns.light_palette("green", n_colors=Labels.unique().size)
                # lut = dict(zip(map(str, Labels.unique()), pal))
                pal = palL
                lut = lutL
                cList.append(pd.Series(Labels.astype('str'), index=daPlotpd.index).map(lut))  # Mapping colours to rows
                legendList.append((Labels.unique().sort_values(), lut))

            elif dim in mDims:
#                 pal = sns.diverging_palette(220, 20, n=Labels.unique().size)
#                 lut = dict(zip(map(str, Labels.unique().sort_values()), pal))  # Note .sort_values() required here.
                pal = palM
                lut = lutM
                cList.append(pd.Series(Labels.astype('str'), index=daPlotpd.index).map(lut))  # Mapping colours to rows
                legendList.append((Labels.unique().sort_values(), lut))

            # If mapping symmetries, use flattened list set above for colour mapping
            # This will keep things consistent over all sym sub-labels
            # NOTE - symFlag setting allows for case when dimMap['Sym'] is missing (will throw an error otherwise)
            # 24/06/20 - removed "and" here - throws error for unmapped types. Nope, reinstated - causes other issues!
            # 14/01/21 - hopefully fixed by using symList, and setting empty default.
            elif symFlag and (dim in dimMap['Sym']):
            # elif symFlag and (dim in symList):
            # elif symFlag or (dim in symList):
#                 pal = sns.color_palette("hls", np.unique(symList).size)
#                 lut = dict(zip(map(str, np.unique(symList)), pal))

                pal = palSym
                lut = lutSym
                cList.append(pd.Series(Labels.astype('str'), index=daPlotpd.index).map(lut))  # Mapping colours to rows
                legendList.append((Labels.unique().sort_values(), lut))
                # pass

            # For other dims, set to categorical mapping.
            # Should check for shared dims here, e.g. Syms, for mapping...?
            else:
#                 pass
#                 print(dim)
#                 print(Labels.unique().size)
#                 pal = sns.color_palette("hls", Labels.unique().size)
                pal = sns.color_palette("dark", Labels.unique().size)
#                 pal = sns.light_palette("blue", n_colors=Labels.unique().size)
                lut = dict(zip(map(str, Labels.unique()), pal))
                cList.append(pd.Series(Labels.astype('str'), index=daPlotpd.index).map(lut))  # Mapping colours to rows
                legendList.append((Labels.unique().sort_values(), lut))

            # For legend, keep also Labels.unique, lut, pairs
#             legendList.append((Labels.unique(), lut, dim))  # Include dim... redundant if Labels is pd.index with name
#             legendList.append((Labels.unique(), lut))


        # Stack to dataframe
        colors = pd.DataFrame(cList[0])
        for c in cList[1:]:
            colors = colors.join(pd.DataFrame(c))

        # Rounding for long labels on xDim
        # if daPlotpd.columns.name in ['Eke','Ehv']:
        # NOW SET IN multiDimXrToPD()

        # *** Plot with modified clustermap code.
        g = snsMatMod.clustermap(daPlotpd,
                  # cmap="vlag", center = 0,
                  # Turn off the clustering
                  row_cluster=False, col_cluster=False,
                  # Add colored class labels
                  row_colors=colors, col_colors=None, # )  # ,
                  # Make the plot look better when many rows/cols
                  # NOTE on ticklabels - pass true/false, auto, int or list for more control.
                  linewidths=0, xticklabels="auto", yticklabels=False,
                  # Some other additional, optional, args...
                  figsize = figsize, cmap = cmap)


        # Add keys for each label - loop over all sets of variables assigned previously as (labels, lut) pairs and set as (invisible) bar plots
        # Using this method it's only possible to set two sets of legends, split these based on n here.
        # Method from https://stackoverflow.com/questions/27988846/how-to-express-classes-on-the-axis-of-a-heatmap-in-seaborn
        # TO DO: For further control over legend layout may need to add dummy items: https://stackoverflow.com/questions/34212241/control-the-number-of-rows-within-a-legend

        # 14/01/21: testing with preset dims & conditionals
        # TO DO: Should be able to just use lutM list and plot global min-max? Currently may have issues with different sized dims?
        # 03/05/21 - testing label rounding here, for float quantities. This fails either for tuple overwrite or on lookup later, ugh. HORRIBLE SHIT.
        #           NOW FORCED WITH enumerate(), but it's still pretty UGLY.  Seems OK in testing, working for legend labels only.
        # TODO: fix with Pandas methods? Should be able to fix for ledgend & axis labels?
        # NOTE: this is set for column labels with multiDimXrToPD() above, should also add index rounding there?
        mDimLabel = None
        for n, item in enumerate(legendList):
            # To avoid issues with long floats, round to 3dp
            # TODO: fix this, currently fails for Euler angles case?
            # 26/09/20: added name check here, hopefully this will fix Euler case...TBC...
            # 03/05/21: retested this... failing for all tried variations, due to either tuple immutability and/or lookups later if values changed.  THIS IS HORRIBLE CODE.
            # Also added above for xDim labels.
            if (item[0].dtype == np.float) or (item[0].name in ['Eke','Ehv','t','Euler']):
                labelsR = np.round(item[0],3).astype('str')
            else:
                labelsR = item[0].astype('str')

            # for label in item[0].astype('str'):
            for ln, label in enumerate(item[0].astype('str')):
            # for label in labels:

                # 03/05/21 - testing label rounding here, this fails on lookup later, ugh.
                # if isinstance(label, np.float):
                #     labelR = np.round(label,3).astype('str')
                # else:
                #     labelR = label  #.astype('str')
                labelR = labelsR[ln]

                if verbose:
                    print(f"Legend entry: {label} for {item}")
#                 label = string(label)
#                 if n%2:

                # 1st legend for l,m,mu
#                 if item[0].name in ['l','m','K','Q','Lambda']:
#                     g.ax_col_dendrogram.bar(0, 0, color=item[1][label],label=label, linewidth=0)
#                 elif item[0].name in ['mu','S','mu0']:  # Skip to avoid repetition, assuming m or Q already plotted on same colour scale. TODO: add conditionals here?
#                     pass
# #                 elif symFlag and (item[0].name is 'Cont'):
# #                     g.ax_row_dendrogram.bar(0, 0, color=item[1][label],label=label, linewidth=0)
#                 elif symFlag and (item[0].name in dimMap['Sym']):  # Skip symmetries
#                 # elif symFlag and (item[0].name in symList):  # Skip symmetries
#                     pass

                # 14/01/21: testing with preset dims & conditionals
                if item[0].name in lDims:
                    g.ax_col_dendrogram.bar(0, 0, color=item[1][label],label=label, linewidth=0)
                    # g.ax_col_dendrogram.bar(0, 0, color=item[1][label],label=labelR, linewidth=0)
                elif item[0].name in mDims:
                    if mDimLabel is None:  # Set first found mDim as source of labels - may cause issues in some cases? Should be able to just use lutM list?
                        mDimLabel = item[0].name
                    if item[0].name == mDimLabel:
                        g.ax_col_dendrogram.bar(0, 0, color=item[1][label],label=f"({label})", linewidth=0)
                        # g.ax_col_dendrogram.bar(0, 0, color=item[1][label],label=f"({labelR})", linewidth=0)

                # elif item[0].name in ['mu','S','mu0']:  # Skip to avoid repetition, assuming m or Q already plotted on same colour scale. TODO: add conditionals here?
                #     pass
#                 elif symFlag and (item[0].name is 'Cont'):
#                     g.ax_row_dendrogram.bar(0, 0, color=item[1][label],label=label, linewidth=0)
                elif symFlag and (item[0].name in dimMap['Sym']):  # Skip symmetries
                # elif symFlag and (item[0].name in symList):  # Skip symmetries
                    pass

                # 2nd legend for other dims.
                else:
                    # g.ax_row_dendrogram.bar(0, 0, color=item[1][label],label=label, linewidth=0)  # Label with value only
                    # g.ax_row_dendrogram.bar(0, 0, color=item[1][label],label=f"{item[0].name}: {label}", linewidth=0)  # Label with category too.
                    g.ax_row_dendrogram.bar(0, 0, color=item[1][label],label=f"{item[0].name}: {labelR}", linewidth=0)  # Label with category too.

        # Add symmetry labels separately from master list to avoid repetition
        if symFlag:
            for label in np.unique(symList):
                g.ax_row_dendrogram.bar(0, 0, color=lutSym[label],label=label, linewidth=0)

#             g.ax_col_dendrogram.legend(title=n, loc="center") #, ncol=n+1) # , bbox_to_anchor=(0.47, 0.8), bbox_transform=plt.gcf().transFigure)

        # Add legends for the bar plots
        ncol = 2  #np.unique(daPlot.LM.l).size
        if 'ADM' in daPlot.dims:
            QNlabels = 'K,(Q,S)'
        else:
            QNlabels = 'l,(m,mu)'

        g.ax_col_dendrogram.legend(title=QNlabels, loc="center", ncol = ncol, bbox_to_anchor=(0.1, 0.6), bbox_transform=plt.gcf().transFigure)
        g.ax_row_dendrogram.legend(title='Categories', loc="center", ncol = 2, bbox_to_anchor=(0.1, 0.4), bbox_transform=plt.gcf().transFigure)

        # Additional anootations etc.
        # Plot titles: https://stackoverflow.com/questions/49254337/how-do-i-add-a-title-to-a-seaborn-clustermap
        # g.fig.suptitle(f"{daPlot.attrs['file']}, pType={pType}, thres={thres}")  # Full figure title
        if thres is not None:
            thresStr = np.round(thres, 2)
        else:
            thresStr = 'None'

        # Set main plot title if not passed.
        if titleString is None:
            if hasattr(daPlot,'jobLabel'):
                titleString = daPlot.attrs['jobLabel']
            else:
                titleString = daPlot.attrs['file']

        g.ax_heatmap.set_title(f"{titleString}\n Plot type = {daPlot.attrs['pTypeDetails']['Type']}, threshold = {thresStr}\n Inc. cross-section {SFflag}, log10 {logFlag}")  # Title heatmap subplot

        # sns.reset_orig()  # Reset gloabl plot params - leaves rc settings, also triggers Matplotlib errors
        sns.reset_defaults()

        return daPlot, daPlotpd, legendList, g

    return daPlot
