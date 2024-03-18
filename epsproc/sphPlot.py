# -*- coding: utf-8 -*-
"""
ePSproc spherical polar plotting functions

Collection of functions for plotting 2D spherical polar data :math:`I(\\theta, \\phi)`.

Main plotters are:

* :py:func:`sphSumPlotX` for arrays of :math:`I(\\theta,\\phi)`
* :py:func:`sphFromBLMPlot` for arrays of :math:`\\beta_{LM}` parameters.

28/11/19
            * Changed plotTypeSelector to dictionary method (and added more methods).

13/09/19

            * Added sphFromBLMPlot().
            * Some additional plot type fuctionality added.
            * Some tidying up and reorganising... hopefully nothing broken...

14/08/19    v1


TODO
----

- Code for switching backend plotter. (Rudiments in place 13/09/19.)
- Generalise plotting for more dimensions.
- More sophisticated/flexible Matplotlib implementation.
- Get Holoviews plotting working.

"""

# imports
import numpy as np
import xarray as xr

# For plotting functions
import matplotlib.pyplot as plt
from matplotlib import cm, colors

# Patch for 3D plotter changes for MPL > v3.4
import matplotlib
if float(matplotlib.__version__[:-2]) < 3.4:
    from mpl_toolkits.mplot3d import Axes3D
else:
    Axes3D = 'Not required'

import copy    # For attrs deepcopy.

# Plotly (optional)
try:
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots
except ImportError as e:
    if e.msg != "No module named 'plotly'":
        raise
    print('* plotly not found, plotly plots not available. ')

# Local functions
from epsproc.sphCalc import sphCalc
from epsproc.util.env import isnotebook
from epsproc.util.misc import checkDims

#***** Plotting top-level functions & logic

# Define plot types
def plotTypeSelector(dataPlot = None, pType = 'a', axisUW = 'Eke', returnDict = False):
    """
    Set plotting data type.

    Parameters
    ----------
    dataPlot : np.array, Xarray
        Data for plotting, output converted to type defined by pType.
        Note: this parameter is required, unless returnDict = True.

    pType : char, optional, default 'a'
        Set type of plot.

        * 'a' (abs)       = np.abs(dataPlot)
        * 'a2' (abs^2)    = np.abs(dataPlot**2)
        * 'r' (real)      = np.real(dataPlot)
        * 'i' (imag)      = np.imag(dataPlot)
        * 'p' (product)   = dataPlot * np.conj(dataPlot)
        * 'phase'         = np.angle(dataPlot)
        * 'phaseUW'       = np.unwrap(np.angle(dataPlot), axis = axisUW)

    axisUW : str, optional, default 'Eke'
        Axis to use for phase unwrapping (for pType = 'phaseUW').
        Axis name must be in passed Xarray.

    returnDict : optional, default = False
        If true, return dictionary of types & methods instead of data array.

    Returns
    -------

    dataPlot : Xarray or np.array
        Input data structure converted to pType.

    pTypeDict : dictionary
        Structure of valid types, formula and functions.

    """

    # # Set dataPlot.values, otherwise attrs are dropped.
    # if pType == 'a':
    #     dataPlot.values = np.abs(dataPlot)
    # elif pType == 'r':
    #     dataPlot.values = np.real(dataPlot)
    # elif pType == 'i':
    #     dataPlot.values = np.imag(dataPlot)
    # elif pType == 'p':
    #     # TODO: check conj defns with Xarray
    #     dataPlot.values = np.abs(dataPlot * np.conj(dataPlot))
    # elif pType == 'a2':
    #     dataPlot.values = np.abs(dataPlot**2)
    # elif pType == 'phase':
    #     dataPlot.values = np.angle(dataPlot)
    # elif pType == 'phaseUW':
    #     axisNum = dataPlot.get_axis_num(axisUW)
    #     dataPlot.values = np.unwrap(np.angle(dataPlot), axis = axisNum)

    # 28/11/19 Set plot type as dictionary
    # Use this instead of if/el list above, and can also be returned for further use.
    # Implement fns as lambdas for flexibility
    pTypeDict = {'a':   {'Type':'abs',  'Formula':'np.abs(dataPlot)',   'Exec': lambda x: np.abs(x)},
                 'a2':  {'Type':'abs(^2)','Formula':'np.abs(dataPlot**2)','Exec': lambda x: np.abs(x**2)},
                 'a2b': {'Type':'(abs)^2','Formula':'np.abs(dataPlot)**2','Exec': lambda x: np.abs(x)**2},
                 'r':   {'Type':'real', 'Formula':'np.real(dataPlot)',  'Exec': lambda x: np.real(x)},
                 'i':   {'Type':'imaginary', 'Formula':'np.imag(dataPlot)',  'Exec': lambda x: np.imag(x)},
                 'p':   {'Type':'product (x*conj(x))', 'Formula':'dataPlot * np.conj(dataPlot)',  'Exec': lambda x: x * np.conj(x)},
                 'pa':   {'Type':'abs(x*conj(x))', 'Formula':'dataPlot * np.conj(dataPlot)',  'Exec': lambda x: np.abs(x * np.conj(x))},
                 'pr':   {'Type':'real(x*conj(x))', 'Formula':'dataPlot * np.conj(dataPlot)',  'Exec': lambda x: np.real(x * np.conj(x))},
                 'pi':   {'Type':'imag(x*conj(x))', 'Formula':'dataPlot * np.conj(dataPlot)',  'Exec': lambda x: np.imag(x * np.conj(x))},
                 'phase':{'Type':'phase', 'Formula':'np.angle(dataPlot)',  'Exec': lambda x: np.angle(x)},
                 'phaseUW':{'Type':'phase (unwrapped)', 'Formula':'np.unwrap(np.angle(dataPlot), axis = axisUW)',  'Exec': lambda x: np.unwrap(np.angle(x), axis = axisNum)}
                 }

    if returnDict:
        return pTypeDict

    # Set pType in output plotting data Xarray
    if pType == 'phaseUW':
        axisNum = dataPlot.get_axis_num(axisUW)  # Special case with axis set

    dataPlot.values = pTypeDict[pType]['Exec'](dataPlot)

    dataPlot.attrs['pType'] = pType
    dataPlot.attrs['pTypeDetails'] = pTypeDict[pType]

    return dataPlot

# Plot a set of mfpads using Holoviews
# TODO: finish this!!!
def sphPlotHV(dataIn):
    # Remove multilevel indexes, these casuse issues with HV
    dataPlot = dataIn.copy()
    if 'Sym' in dataPlot.dims:      # This is a bit dodgy, and may fail in some cases.
        dataPlot['Sym'] = dataPlot['Cont']

    if 'Euler' in dataPlot.dims:
        dataPlot['Euler'] =  np.arange(0,dataPlot.Euler.size)

    # Convert sph to cart coords
    # ACTUALLY, this will be a bit tricky at da level - best to set a new da for (x,y,z) outputs...?
    # dataPlot = dataPlot * dataPlot.conj()
    # theta, phi = np.meshgrid(dataPlot.Theta, dataPlot.Phi)
    # X = dataPlot * np.sin(theta) * np.cos(phi)
    # Y = dataPlot * np.sin(phi) * np.sin(theta)
    # Z = dataPlot * np.cos(phi)




# Plot MFPADs from a set of BLM
def sphFromBLMPlot(BLMXin, res = 50, pType = 'a', plotFlag = False, facetDim = None, backend = 'mpl', fnType = None, conj = False, **kwargs):
    r'''
    Calculate spherical harmonic expansions from BLM parameters and plot.

    Surfaces calculated as:

    .. math::
        I(\theta,\phi)=\sum_{L,M}\beta_{L,M}Y_{L,M}(\theta,\phi)


    Parameters
    ----------
    dataIn : Xarray
        Input set of BLM parameters, or other (L,M) exapansion dataType.

    res : int, optional, default 50
        Resolution for output (theta,phi) grids.

    pType : char, optional, default 'r' (real part)
        Set (data) type to plot. See :py:func:`plotTypeSelector`.

    plotFlag : bool, optional, default False
        Set plotting True/False.  Note that this will plot for all facetDim.

    facetDim : str, optional, default None
        Dimension to use for subplots.
        Currently set for a single dimension only.
        For matplotlib backend: one figure per surface.
        For plotly backend: subplots per surface.

    backend : str, optional, default 'mpl' (matplotlib)
        Set backend used for plotting. See :py:func:`sphSumPlotX` for details.

    fnType : str, optional, default = 'sph'
        Set function for expansion parameters, default is YLM from scipy.special.sph_harm.
        See :py:func:`ep.sphCalc` for details.

    conj : bool, optional, default = False
        Use conjugate harmonics.

    **kwargs
        Additional args passed to backend plotting routines.

    Returns
    -------
    Xarray
        Xarray containing the calcualted surfaces I(theta,phi,...)

    fig
        List of figure handles.

    '''

    #*** Check dataType
    BLMX = BLMXin.copy()  # Default case

    if 'dataType' in BLMXin.attrs.keys():
        # For ADMs, can use this plotting routine if S=0 only.
        # NOTE - this should work if other singleton dims are present, but will need to explicitly set facetDims for plots, otherwise will throw an error
        if (BLMXin.attrs['dataType'] == 'ADM'):
            if all(BLMXin.S == 0):
                # NOTE: Squeeze here to kill S dim, but this will also drop (l,m) if singleton. Needs a fix.
                # TODO: fix this!
                # UPDATE 05/01/22: in testing with XR v0.15 this is OK without squeeze(). Have left in with additional dim check for now/safety.
                #                   EDIT: actually better with .squeeze('S') to be explict and ensure OK for arb dims.
                BLMX = BLMXin.copy().unstack('ADM').rename({'K':'l', 'Q':'m'}).drop('S').stack({'BLM':('l','m')}).squeeze('S')

                # if BLMX.l.size > 1:
                #     BLMX = BLMX.squeeze()

            else:
                print('***ADM set with non-zero S, skipping Ylm routine.')
                BLMX = None

        # For normal BLMs, no action required.
        elif BLMXin.attrs['dataType'] == 'BLM':
            # BLMX = BLMXin
            pass

        # Default case - set as BLM case but include dataType warning.
        else:
            print(f"*** sphPlot dataType = {BLMXin.attrs['dataType']} not recognised, trying anyway.")
            # BLMX = BLMXin

    else:
        print("*** sphPlot no dataType set. Trying anyway. Set data.attrs['dataType'] for more control.")
        # BLMX = BLMXin

    # Ensure attrs are deepcopied
    BLMX.attrs = copy.deepcopy(BLMXin.attrs)  # THIS WORKS ALSO FOR NESTED DICT CASE
                                            # Definitely required in XR2022.3, 2022.6, but should be fixed in later versions, see https://github.com/pydata/xarray/issues/2835


    fig = None
    if BLMX is not None:
        # Calculate YLMs
        if hasattr(BLMX,'normType') and (fnType is None):
            fnType = BLMX.attrs['normType']
            print(f'Using {fnType} betas (from BLMX array).')
        elif hasattr(BLMX,'harmonics') and (fnType is None):
            fnType = BLMX.attrs['harmonics']['kind']
            print(f'Using {fnType} betas (from BLMX array).')
        elif fnType is not None:
            print(f'Using {fnType} betas (as passed).')
        else:
            fnType = 'sph'
            print(f'Using default {fnType} betas.')


        YLMX = sphCalc(BLMX.l.max(), res=res, fnType=fnType, conj = conj)

        if 'BLM' in BLMX.dims:
            YLMX = YLMX.rename({'LM':'BLM'})    # Switch naming for multiplication & plotting

        # Calculate MFPADs (theta,phi)
        dataPlot = BLMX*YLMX

        # 27/09/22 VERY basic dim name handling added here, but should generalise.
        if 'BLM' in dataPlot.dims:
            dataPlot = dataPlot.rename({'BLM':'LM'})    # Switch naming back for plotting function

        # dataPlot.attrs = BLMX.attrs # Ensure attrs propagated
        # dataPlot.attrs.update(YLMX.attrs)  # Add YLMX attrs - may overwrite
        dataPlot.attrs = YLMX.attrs # Ensure attrs propagated
        dataPlot.attrs.update(BLMX.attrs)  # Add BLMX attrs - may overwrite
        if 'harmonics' in dataPlot.attrs.keys():
            dataPlot.attrs['harmonics'].update(YLMX.attrs['harmonics'])  # Ensure YLM settings propagated correctly.

        dataPlot.attrs['dataType'] = 'Itp'
        dataPlot.attrs['long_name'] = r'I(\theta,\phi)'
        dataPlot.attrs['normType'] = fnType

        # Pass data to plotting function
        if plotFlag:
            fig = sphSumPlotX(dataPlot, pType = pType, facetDim = facetDim, backend = backend, **kwargs)
        else:
            fig = None

        # 27/03/23 - changed to always run plotter, allows for figure return.
        # 12/07/23 - removed and modified conditional case above instead, otherwise can get issue with some cases for calc only where not all args passed.
        # fig = sphSumPlotX(dataPlot, pType = pType, facetDim = facetDim, backend = backend, plotFlag = plotFlag, **kwargs)


    return dataPlot, fig


# Sum and plot spherical functions from an Xarray.
# TODO: This currently assumes matplotlib cm is already loaded.
# TODO: More plot types.
# TODO: More careful testing - not totally sure if summation & axes consistent here.
# TODO: add dim checks, see classes._plotters.padPlot for some existing implementations.
def sphSumPlotX(dataIn, pType = 'a', facetDim = 'Eke', surfMap = 'R', backend = 'mpl',
                convention = 'phys', titleString = None, plotFlag = True, verbose = True, axisUW = None,
                **kwargs):
    '''
    Plot sum of spherical harmonics from an Xarray.

    Parameters
    ----------
    dataIn : Xarray
        Input structure can be

        * Set of precalculated Ylms, dims (theta,phi) or (theta,phi,LM).
        * Set of precalculated mfpads, dims (theta,phi), (theta,phi,LM) or (theta,phi,LM,facetDim).
        * If (LM) dimension is present, it is summed over before plotting.
        * If facetDim is present this is used for subplots, currently only one facetDim is supported here.

    pType : char, optional, default 'a' (abs value)
        Set (data) type of plot. See :py:func:`plotTypeSelector`.

    facetDim : str, optional, default Eke
        Dimension to use for subplots.

        * Currently set for a single dimension only.
        * For matplotlib backend: one figure per surface.
        * For plotly backend: subplots per surface.

    surfMap : str, optional, default = None
        Additional specification to use for surface colour specification.
        NOTE: currently only implemented for Plotly backend (via surfacecolor).
        If None not specified == z surf map.
        If R == radial value surf map.

    backend : str, optional, default 'mpl'
        Set backend used for plotting.

        - mpl   matplotlib: basic 3D plotting, one figure per surface.
        - pl    plotly: fancier 3D plotting, interactive in Jupyter but may fail at console.
                Subplots for surfaces.
        - hv    holoviews: fancier plotting with additional back-end options.
                Can facet on specific data types.

    convention : str, optional, default = 'phys'
        Spherical polar coord convention, see :py:func:`epsproc.sphToCart`

    titleString : str, optional, default = None
        Additional info to use for plot title.

    plotFlag : bool, optional, default = True
        Set plotFlag=False bypass for plotter object return only. (Required for Plotly backend only.)

    axisUW : str, optional, default = None
        Unwrap axis for phase maps, only used if pType = 'phaseUW'

    **kwargs
        Additional args passed to backend plotting routines.

    Returns
    -------
    fig
        List of figure handles.

    Examples
    --------
    >>> YlmX = sphCalc(2, res = 50)
    >>> sphSumPlotX(YlmX)

    Note
    ----
    Pretty basic functionality here, should add more colour mapping options and multiple plots, alternative back-ends, support for more dimensions etc.

    More sophisticated dim handling is implemented in padPlot(), but for class only.


    '''

    # Sum over QNs if necessary
    # if (dataIn.LM.shape[0]) > 1:
    if 'LM' in dataIn.dims:
        dataPlot = dataIn.sum(dim='LM')
    else:
        dataPlot = dataIn

    # (theta,phi) grid from Xarray coords
    theta, phi = np.meshgrid(dataPlot.Theta, dataPlot.Phi)

    if dataPlot.min() < 0:
        print(f"*** WARNING: plot dataset has min value < 0, min = {dataPlot.min().values}. This may be unphysical and/or result in plotting issues.")

    # Set data according to type of plot selected
    dataPlot = plotTypeSelector(dataPlot, pType = pType, axisUW = axisUW)

    # Set main plot title if not passed.
    if titleString is None:
        if hasattr(dataIn,'jobLabel'):
            titleString = dataIn.attrs['jobLabel']
        elif hasattr(dataIn,'file'):
            titleString = dataIn.attrs['file']
        else:
            titleString = ""

    if verbose:
        print(f"Sph plots: {titleString}")  # Throw out titlestring here for ease.

        # Switch plot function based on backend
        # print('Plotting with {0}'.format(backend))
        print(f"Plotting with facetDims={facetDim}, pType={pType} with backend={backend}.")

    fig = []

    # Matplotlib
    if backend == 'mpl':
        # Check dimensionality - loop over facetDim if necessary
        # Return list of fig handles
        if len(dataPlot.dims) > 2:

            if verbose:
                print('Data dims: {0}, subplots on {1}'.format(dataPlot.dims, facetDim))

            # Basic loop over facetDim
            # Fails for dims with multi-index, or repeated values for sub-indexes if used as selector.
            # Therefore, use positional numerical index...
            # for nPlot in dataPlot[facetDim]:
            #     fig.append(sphPlotMPL(dataPlot.sel({facetDim:[nPlot]}).squeeze(), theta, phi))  # Will fail if dims>2 passed.

            # Loop over facetDim with full dataArray as selector, allows for MultiIndex cases.
            for nPlot in dataPlot[facetDim]:
                # fig.append(sphPlotMPL(dataPlot.sel({facetDim:nPlot}).squeeze(), theta, phi,  convention = convention, tString = f"{facetDim}: {nPlot.item()}"))
                # {dataPlot[facetDim][nPlot].item()}"))  # Will fail if dims>2 passed.
                fig.append(sphPlotMPL(dataPlot.sel({facetDim:nPlot}).squeeze(), theta, phi,  convention = convention, tString = titleString + f"\n{facetDim}: {nPlot.item()}", **kwargs))

        else:
            # Call matplotlib plotting fn., single surface
            fig.append(sphPlotMPL(dataPlot, theta, phi,  convention = convention, tString = titleString))

    # Plotly - note that faceting is handled directly by Plotly in this case.
    if backend == 'pl':
        fig.append(sphPlotPL(dataPlot, theta, phi, facetDim, surfMap=surfMap, convention = convention, plotFlag = plotFlag, verbose = verbose, **kwargs))

    return fig


#******** Low-level plotting functions

# Set cart coords from spherical polar coords
def sphToCart(R, theta, phi, convention = 'phys', returnType = 'np'):
    r"""
    Convert spherical polar coords :math:`(R,\theta,\phi)` to Cartesian :math:`(X,Y,Z)`.

    Parameters
    ----------
    R, theta, phi : np.arrays
        Spherical polar coords :math:`(R,\theta,\phi)`.

    convention : str, optional, default = 'phys'
        Specify choice of Spherical Polar coordinate system, 'phys' or 'maths' (see note below).

    returnType : str, optional, default = 'np'
        - 'np' return numpy arrays X,Y,Z (default)
        - 'xr' return Xarray Dataset

    Returns
    -------
    X, Y, Z : np.arrays or Xarray Dataset
        Cartesian coords :math:`(X,Y,Z)`.


    Definitions
    ------------
    Conversion defined with the `usual physics convention <https://en.wikipedia.org/wiki/Spherical_coordinate_system#Conventions>`_ , where:

    * :math:`R` is the radial distance from the origin
    * :math:`\theta` is the polar angle (defined relative to the z-axis), :math:`0\leq\theta\leq\pi`
    * :math:`\phi` is the azimuthal angle (defined relative to the x-axis), :math:`0\leq\theta\leq2\pi`

    * :math:`X = R * np.sin(phi) * np.cos(theta)`
    * :math:`Y = R * np.sin(phi) * np.sin(theta)`
    * :math:`Z = R * np.cos(phi)`

    Specify convention = 'maths' to use alternative definition with theta and phi swapped.

    """

    if convention == 'maths':
        # Mathematics convention
        X = R * np.sin(phi) * np.cos(theta)
        Y = R * np.sin(phi) * np.sin(theta)
        Z = R * np.cos(phi)

    elif convention == 'phys':
        # Physics convention
        X = R * np.sin(theta) * np.cos(phi)
        Y = R * np.sin(theta) * np.sin(phi)
        Z = R * np.cos(theta)

    else:
        print(f"*** Coordinate convention {convention} not supported.")
        return None

    # if returnType == 'np':

    if returnType == 'xr':
        return xr.Dataset({'X':X,'Y':Y,'Z':Z})
    else:
        return X, Y, Z


# Plot using matplotlib - legacy
# TODO - remove this.
# def sphPlot(dataPlot, theta, phi):
#     sphPlotMatplotLib(dataPlot, theta, phi)

# Plot using matplotlib
def sphPlotMPL(dataPlot, theta, phi, convention = 'phys', tString = None, **kwargs):
    '''
    Plot spherical polar function (R,theta,phi) to a Cartesian grid, using Matplotlib.

    Parameters
    ----------
    dataPlot : np.array or Xarray
        Values to plot, single surface only, with dims (theta,phi).

    theta, phi : np.arrays
        Angles defining spherical polar grid, 2D arrays.

    convention : str, optional, default = 'phys'
        Set spherical polar coord convention, see :py:func:`epsproc.sphToCart`.

    tString : str, optional, default = None
        Text to be used for plot title.
        This will be appended with other data info, if set.
        If facetDim is passed here, this will be used to set the label.

    **kwargs
        Unused, just for calling func convenience.

    Returns
    -------
    fig
        Handle to matplotlib figure.


    '''

    X, Y, Z = sphToCart(dataPlot, theta, phi, convention = convention)

    # Plot in a new figure using Matplotlib
    # 13/07/23 - patched for MPL changes for v3.4+
    # For MPL v3.7+ old code won't work at all, for 3.4-3.6 get a warning, "Axes3D(fig) adding itself to the figure is deprecated since 3.4. Pass the keyword argument auto_add_to_figure=False and use fig.add_axes(ax) to suppress this warning. The default value of auto_add_to_figure will change to False in mpl3.5 and True values will no longer work in 3.6.  This is consistent with other Axes classes."
    if float(matplotlib.__version__[:-2]) < 3.4:
        # Old style
        fig = plt.figure()
        ax = Axes3D(fig)
    else:
        # New style, per https://matplotlib.org/stable/plot_types/3D/surface3d_simple.html#sphx-glr-plot-types-3d-surface3d-simple-py
        fig, ax = plt.subplots(subplot_kw={"projection": "3d"})



    # ax.plot_surface(x, y, z, rstride=1, cstride=1, cmap=cm.jet)
    ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap=cm.jet)
    # ax.axis('equal') # Not implemented for 3D axes
    # Rescale axis to equal, hack from https://stackoverflow.com/questions/8130823/set-matplotlib-3d-plot-aspect-ratio
    scaling = np.array([getattr(ax, 'get_{}lim'.format(dim))() for dim in 'xyz'])
    # ax.auto_scale_xyz(*[[np.min(scaling), np.max(scaling)]]*3)  # Not quite right - can have lopsided axes for asymmetric distributions.
    mRange = np.max([np.abs(np.min(scaling)), np.max(scaling)])
    ax.auto_scale_xyz(*[[-mRange, mRange]]*3)

    # Force axis off
    # ax.set_axis_off()

    #*** Title string
    # Check if facetDim was passed
    if tString is not None:
        if tString in dataPlot.dims:
            facetDim = tString
            tString = f"{facetDim}: {dataPlot['facetDim'].item()}"

        # else:

    elif tString is None:
        # Title with data details.
        if hasattr(dataPlot, 'jobLabel'):
            tString = dataPlot.jobLabel
        elif hasattr(dataPlot, 'file'):
            tString = dataPlot.file
        else:
            tString = 'None test'


    plt.title(tString)
    # plt.show()
    # fig.show()

    return fig


# Plot as Plotly subplots, as function of selected dim.
# Currently set for subplotting over facetDim, but assumes other dims (theta,phi)
def sphPlotPL(dataPlot, theta, phi, facetDim = 'Eke', surfMap = None,
                showbackground = False, showaxes = True,
                rc = None, nCols = 4, norm = 'local', padding = 0.05, camR = 0.85, #'global',
                convention = 'phys', plotFlag = True, verbose = False, **kwargs):
    '''
    Plot spherical polar function (R,theta,phi) to a Cartesian grid, using Plotly.

    Parameters
    ----------
    dataPlot : np.array or Xarray
        Values to plot, single surface only, with dims (theta,phi).

    theta, phi : np.arrays
        Angles defining spherical polar grid, 2D arrays.

    facetDim : str, default 'Eke'
        Dimension to use for faceting (subplots), currently set for single dim only.

    surfMap : str, optional, default = None
        Additional specification to use for Plotly surfacecolor specification.
        If None not specified == z surf map.
        If R == radial value surf map.

    showbackground : bool, default = False
        Show background axis grid and planes if True.

    showaxes : bool, default = True
        Plot axes lines if True.

    rc : array, optional, default = None
        If set, use to define layout grid [rows, columns].
        If not set, this will be set for nCols = 4

    nCols : int, optional, default = 4
        Number of columns in plot layout grid (if required).
        Layout grid will be set as rc = [nData/nCols, nCols]
        Where nData is number of plots.

    norm : str, optional, default = 'global'
        Set for plot normalisation.
        - 'global' : use same (x,y,z) limits for all plots.
        - 'local' : auto set (x,y,z) limits for each plot.

    padding : float, optional, default = 0.05
        Axis padding for plot scaling, as %age of R max.

    camR : float, optional, default = 0.85
        Use to move camera to zoom in on plots, with camera = dict(eye=dict(x=camR, y=camR, z=camR)).
        Plotly default case (1,1,1) usually doesn't fill render box.

    convention : str, optional, default = 'phys'
        Spherical polar coord convention, see :py:func:`epsproc.sphToCart`

    plotFlag : bool, optional, default = True
        Set plotFlag=False bypass for plotter object return only.

    **kwargs
        Unused, just for calling func convenience.


    Returns
    -------
    fig
        Handle to figure.

    Notes
    -----
    For additional dev notes, see:
    - http://localhost:8889/notebooks/dev/ePSproc/plottingDev/plotly_subplot_tests_051020.ipynb
    - http://localhost:8889/notebooks/dev/ePSproc/classDev/ePSproc_multijob_class_tests_N2O_011020_Stimpy_PLrenderTESTING_local.ipynb

    For Jupyter use: currently (Oct. 2020) working only for Notebook Export to HTML (not nbsphinx), and max of 12 subplots. Possible issues with rendering in Firefox (80.0.1, Oct. 2020).
    Should be fixable with some effort/testing, see https://plotly.com/python/renderers/

    TODO:

    - More playing around with Plotly.
    - Camera control and linking, e.g. https://community.plotly.com/t/synchronize-camera-across-3d-subplots/22236
    - Update 14/03/23: added axis and background plotting options. Need to consolidate with global styles.
    - Update 16/03/23: better facet dim handling, and basic camera control added.

    For JupyterLab, need additional extensions - see https://plotly.com/python/getting-started/#jupyterlab-support:
    - `conda install -c conda-forge -c plotly jupyter-dash`
    - `jupyter labextension install jupyterlab-plotly`

    In some cases may get partially working installation with, e.g., blank surface plots, or plots via HV only. This usually means JupyterLab needs a restart (and maybe a rebuild).
    For more see https://plotly.com/python/troubleshooting/

    '''

    # Set up subplots
    # 31/03/22 type check to fix issues with padPlot wrapper and None type facetDims.
    # 12/08/22 also fixed plot titles & X,Y,Z clac for None case below.
    if facetDim is not None:
        nData = dataPlot[facetDim].size
    else:
        # facetDim = 'Eke'
        # nData = dataPlot[facetDim].size
        nData = 1
        norm = 'global'  # Force to global, otherwise norm-per-facet will fail.

        if rc is None:
            rc = [1,1]

    if rc is None:
        # nCols = 6
        if nData < nCols:   # Set for single row layout
            nCols = nData
        rc = [nData/nCols, nCols]

    # rc = np.round(rc).astype(np.int)  # Fix dtype - Plotly throws type error for np types however.
    rc = [int(np.ceil(rc[0])), int(np.ceil(rc[1]))]

    pType = {'type':'surface'}
    specs = [[pType] * rc[1] for i in range(rc[0])]  # Set specs as 2D list of dicts.

    # Check/set some global props - assumes Xarray datatype at the moment (also assumed later!)
    # 16/03/23: updated sphToCart for Dataset return, and this code to use Dataset with coords.
    if isinstance(dataPlot, xr.core.dataarray.DataArray):
        # Check max coord limit, to use for all subplot ranges
        XYZds =sphToCart(dataPlot,dataPlot.Theta,dataPlot.Phi, returnType='xr')  # NOTE - this is NOT used for plotting, just global limits.

        # Get max & mins per coord
        XYZds.attrs['dimMax'] = XYZds.max().to_array()
        XYZds.attrs['dimMin'] = XYZds.min().to_array()

        # Max & min global
        XYZds.attrs['CMax'] = XYZds.attrs['dimMax'].max()
        XYZds.attrs['CMin'] = XYZds.attrs['dimMin'].min()
        XYZds.attrs['Rmax'] = np.max([np.abs(XYZds.attrs['CMin']), XYZds.attrs['CMax']])

        # Axes settings global case
        # padding = 0.1 * XYZds.attrs['Rmax']  # Padding %age
        XYZds.attrs['aRanges'] = dict(range=[-(XYZds.attrs['Rmax'] + padding*XYZds.attrs['Rmax']),
                                        XYZds.attrs['Rmax']+padding*XYZds.attrs['Rmax']])

        # Set facet properties if required
        # Set subplot titles, see https://stackoverflow.com/questions/53991810/how-to-add-subplot-title-to-3d-subplot-in-plotly
        if facetDim is not None:
            # Check dims are OK
            dimDict = checkDims(XYZds, refDims=facetDim)

            # Per subplot settings... max/min per facet dim.
            XYZds.attrs['facetMax'] = XYZds.max(dim=dimDict['extra'])
            XYZds.attrs['facetMin'] = XYZds.min(dim=dimDict['extra'])
            XYZds.attrs['CfacetMax'] = XYZds.attrs['facetMax'].to_array().max('variable')
            XYZds.attrs['CfacetMin'] = XYZds.attrs['facetMin'].to_array().min('variable')
            XYZds.attrs['RfacetMax'] = xr.concat([XYZds.attrs['CfacetMax'], np.abs(XYZds.attrs['CfacetMin'])], dim='m').max('m')
            XYZds.attrs['RfacetMaxRel'] = XYZds.attrs['RfacetMax'].max()/XYZds.attrs['RfacetMax']

            # Axes settings facet case
            # padding = 0.1 * XYZds.attrs['Rmax']  # Padding %age
            # XYZds.attrs['aRanges'] = dict(range=[-(XYZds.attrs['Rmax'] + padding), XYZds.attrs['Rmax']+padding])
            # Max/min R per facet dim with padding
            XYZds.attrs['aRangesFacet'] = xr.Dataset({'min':-(XYZds.attrs['RfacetMax'] + padding*XYZds.attrs['RfacetMax']),
                                                    'max':XYZds.attrs['RfacetMax']+padding*XYZds.attrs['RfacetMax']}).to_array(dim='m')

            if norm != 'global':
                # Add titles including scale factor
                # Note <br> for newline in Plotly text object
                titles = [f"{facetDim}: {item.item()}<br>Relative scale={np.round(XYZds.attrs['RfacetMaxRel'][n].data,1)}" for n,item in enumerate(dataPlot[facetDim])]
            else:
                titles = [f"{facetDim}: {item.item()}" for item in dataPlot[facetDim]]

        else:
            # titles = [f"{facetDim}: {item.item()}" for item in dataPlot[facetDim]]
            titles = ['']  # Singleton case.

        # # Check max coord limit, to use for all subplot ranges
        # # Should rewrite sphToCart for Xarray output?
        # X,Y,Z =sphToCart(dataPlot,dataPlot.Theta,dataPlot.Phi)  # NOTE - this is NOT used for plotting, just global limits.
        # # X,Y,Z =sphToCart(dataPlot,theta,phi)  # NOTE - (t,p) here predefined 2D arrays!
        #
        # Cmax = np.max([X.max(), Y.max(), Z.max()])  # Check max & min coord values - may be unequal in some cases.
        # Cmin = np.min([X.min(), Y.min(), Z.min()])
        # Rmax = np.max([np.abs(Cmin), Cmax])
        # padding = 0.1 * Rmax  # Padding %age
        # aRanges = dict(range=[-(Rmax + padding), Rmax+padding])
        #
        # # Set subplot titles, see https://stackoverflow.com/questions/53991810/how-to-add-subplot-title-to-3d-subplot-in-plotly
        # if facetDim is not None:
        #     # Max per subplot
        #     dimsPlot = set(dataPlot.dims) - set([facetDim])
        #     Xfacet = X.max(dim=['Theta','Phi'])
        #     Yfacet = Y.max(dim=['Theta','Phi'])
        #     Zfacet = Z.max(dim=['Theta','Phi'])
        #     RmaxFacet = np.max(np.array([Xfacet,Yfacet,Zfacet]), axis=0)
        #     RmaxFacetRel = RmaxFacet.max()/RmaxFacet  # Relative scale factor compared to largest magnitude case.
        #
        #     if norm != 'global':
        #         # Add titles including scale factor
        #         # Note <br> for newline in Plotly text object
        #         titles = [f"{facetDim}: {item.item()}<br>Relative scale={np.round(RmaxFacetRel[n],1)}" for n,item in enumerate(dataPlot[facetDim])]
        #     else:
        #         titles = [f"{facetDim}: {item.item()}" for item in dataPlot[facetDim]]
        #
        # else:
        #     # titles = [f"{facetDim}: {item.item()}" for item in dataPlot[facetDim]]
        #     titles = ['']  # Singleton case.

    # Set up subplots
    fig = make_subplots(rows=rc[0], cols=rc[1], specs=specs, subplot_titles=titles)
    nPlots = rc[0] * rc[1]

    # Add surfaces
    # From https://plot.ly/python/3d-subplots/

    n=0
    for rInd in range(1,rc[0]+1):
        for cInd in range(1,rc[1]+1):

            if verbose:
                print(f"*** Plotting for [{rInd},{cInd},{n}]")

            if n < nData:

                # Set data (NOTE - now repeats above in Xarray case)
                if facetDim is not None:
                    R = dataPlot.sel({facetDim:dataPlot[facetDim][n]})
                    # X,Y,Z = sphToCart(dataPlot.sel({facetDim:dataPlot[facetDim][n]}),theta,phi)  # Bit ugly, probably a better way to select here.
                else:
                    R = dataPlot
                    # X,Y,Z = sphToCart(dataPlot,theta,phi)  # Singleton case

                X,Y,Z = sphToCart(R,theta,phi)

                # Additional surface colour map options
                if surfMap is None:
                    surfMapArray = Z
                elif surfMap == 'R':
                    surfMapArray = R

                # Set plot object
                showscale = False
                # if rInd == (rc[0]):
                #     showscale = True  # Set scale bar for row?  Only for global scaling? With this code get multiple cbars?

                fig.add_trace(
                    go.Surface(x=X, y=Y, z=Z, colorscale='Viridis', showscale=showscale, surfacecolor=surfMapArray),
                    row=rInd, col=cInd)  # Needs some work here...
                    # title=[facetDim + '=' + str(dataPlot[facetDim][n].data.item())]),

                # fig['layout'].update(scene=dict(aspectmode="data"))  # Try and fix aspect ratio issues - getting stuck on first subplot?  Seems independent of setting here.

                # DOESN'T WORK
                # From https://github.com/plotly/plotly.py/issues/70
                # fig['layout'].update(scene=dict(
                #                     aspectmode='manual',
                #                     aspectratio=go.layout.scene.Aspectratio(
                #                     x=x.ptp(), y=y.ptp(), z=z.pyp())))

                n=n+1

                # Set additional axis properties if required
                # Needs some work, see
                # - https://plotly.com/python/3d-axes/
                # - https://plotly.com/python/3d-surface-plots/
                # - https://plotly.com/python/3d-subplots/#3d-surface-subplots
                # TODO: try linking cameras...? https://community.plotly.com/t/synchronize-camera-across-3d-subplots/22236

                # Set string for "scene" (axis) object to update - will be labelled scene1, scene2... by Plotly.
                sceneN = f'scene{n}'
                if norm == 'global':
                    # Try looping... OK with dict unpacking... huzzah!
                    # NOTE Scene indexing starts at 1, so do this after n increments
                    # options = dict(xaxis = aRanges, yaxis = aRanges, zaxis = aRanges, aspectmode='cube')
                    options = dict(xaxis = XYZds.attrs['aRanges'], yaxis = XYZds.attrs['aRanges'], zaxis = XYZds.attrs['aRanges'], aspectmode='cube')

                else:
                    # options = dict(aspectmode='cube')  # Doesn't seem to work as expected (distorts axes)
                    # options = dict(aspectmode='auto')   # Ah, OK - relative scaling/aspect ratio preserved in this case.
                                                        # BUT zoom set poorly?  Often seems to be x2 smaller than it should?
                                                        # May want to force with per-plot camera options?

                    # With ranges per facet set explicitly
                    aRangesFacet = dict(range=list(XYZds.attrs['aRangesFacet'].sel({facetDim:XYZds.attrs['aRangesFacet'][facetDim][n-1]}).data))
                    options = dict(xaxis = aRangesFacet, yaxis = aRangesFacet, zaxis = aRangesFacet, aspectmode='cube')

                if camR is not None:
                    options['camera'] = dict(eye=dict(x=camR, y=camR, z=camR))

                fig.update_layout(**{sceneN:options})  # No effect of aspect here? auto/cube/data/manual

                # if rInd == (rc[0]):
                #     fig.update_layout(**{sceneN:dict(showscale=True)})  # Set scale bar for row?

                # Remove grid and background?
                if not showbackground:
        #             whiteOptions = dict(backgroundcolor = 'rgb(255,255,255)', color='rgb(255,255,255)')  # White bg
        #             whiteOptions = dict(showbackground=False)   # No bg, keep axes
                    whiteOptions = dict(showbackground=False, visible=False)   # showline=False, showticklabels=False, showaxeslabels=False)   # No bg, no axes
        #             whiteOptions = dict(showbackground=True, showline=False, showticklabels=False)    # zeroline=True)   # No bg, test other options
                    fig.update_scenes(xaxis=whiteOptions, yaxis=whiteOptions, zaxis=whiteOptions)

                # Add (x,y,z) axes only?
                # Use with Rmax or (X,Y,Z) max to define axis extent.
                # Add padding to always show ends of axes.
                if showaxes:
                    if norm == 'global':
                        fig.add_trace(createAxesPL(scale=[XYZds.attrs['Rmax']+padding,XYZds.attrs['Rmax']+padding,XYZds.attrs['Rmax']+padding]),row=rInd, col=cInd)
                                # Max & min global
                                # XYZds.attrs['CMax'] = XYZds.attrs['dimMax'].max()
                                # XYZds.attrs['CMin'] = XYZds.attrs['dimMin'].min()
                                # XYZds.attrs['Rmax'] = np.max([np.abs(XYZds.attrs['CMin']), XYZds.attrs['CMax']])
                                # Axes settings global case
                                # padding = 0.1 * XYZds.attrs['Rmax']  # Padding %age
                                # XYZds.attrs['aRanges'] = dict(range=[-(XYZds.attrs['Rmax'] + padding*XYZds.attrs['Rmax']),
                                #                                 XYZds.attrs['Rmax']+padding*XYZds.attrs['Rmax']])
                    else:
                        # fig.add_trace(createAxesPL(scale=[X.max()+padding,Y.max()+padding,Z.max()+padding]),row=rInd, col=cInd)

                        # TODO: replace above with precalc values from XYZds
                        # Can use XYZds.attrs['aRangesFacet'].sel({facetDim:XYZds.attrs['aRangesFacet'][facetDim][n-1] as above
                        # Or max per (X,Y,Z)
                        # RmaxAxis = XYZds.attrs['RfacetMax'][facetDim][n-1]
                        # fig.add_trace(createAxesPL(scale=),row=rInd, col=cInd)

                        # 16/03/23 - working, but may want to add additional padding to ensure axis visibility?
                        # UPDATE: padding now added. ~5-10% works nicely.
                        axScale = XYZds.attrs['facetMax'].to_array().sel({facetDim:XYZds.attrs['aRangesFacet'][facetDim][n-1]})
                        # print(axScale)
                        fig.add_trace(createAxesPL(scale=axScale+axScale*padding),row=rInd, col=cInd)

    # fig.show()
    # Check if notebook and output
    # Allow for plotFlag=False bypass for object return only.
    if isnotebook() and plotFlag:
        display(fig) # Use IPython display(), works nicely for notebook output
        # Q: is display() always loaded automatically? I think so.
    elif plotFlag:
        fig.show()  # Try fig.show(), although may not work in all cases.

    return fig


def createAxesPL(scale = [1,1,1], color='rgb(0,0,0)', width=4):
    """
    Create (x,y,z) axes lines as Plotly go.Scatter3d object
    Based on https://stackoverflow.com/questions/42301481/adding-specific-lines-to-a-plotly-scatter3d-plot

    Basically set and draw point pairs, then lines based on these

    Parameters
    ----------
    scale : list or array, default = [1,1,1]
        [x,y,z] limits for the axes. Note +/- extents are drawn.
        TODO: added +/- options and dictionary option?

    color : str, default = 'rgb(0,0,0)'
        Define line colour. See Plotly docs for supported values.

    width : int, default = 4
        Line width


    Returns
    -------

    axes : plotly go.Scatter3d object.


    Examples
    --------

    >>> # Draw unit axes
    >>> import plotly.graph_objs as go
    >>> fig = go.Figure(data=createAxesPL())
    >>> fig.show()

    >>> # Draw axes with custom scale
    >>> import plotly.graph_objs as go
    >>> fig = go.Figure(data=createAxesPL([5,10,1]))
    >>> fig.show()

    """

    # Define points
    x = [-scale[0],scale[0],0,0,0,0]
    y = [0,0,-scale[1],scale[1],0,0]
    z = [0,0,0,0,-scale[2],scale[2]]

    # Define pairs (the start and end point for each line)
    pairs = [(0,1), (2,3), (4,5)]

    # Define lines
    x_lines = list()
    y_lines = list()
    z_lines = list()

    # create the coordinate list for the lines
    for p in pairs:
        for i in range(2):
            x_lines.append(x[p[i]])
            y_lines.append(y[p[i]])
            z_lines.append(z[p[i]])
        x_lines.append(None)
        y_lines.append(None)
        z_lines.append(None)

    # Draw lines with Scatter3d
    axes = go.Scatter3d(
        x=x_lines,
        y=y_lines,
        z=z_lines,
        mode='lines',
        line=dict(color=color, width=width),showlegend=False)

    return axes
