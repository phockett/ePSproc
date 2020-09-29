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
# For plotting functions
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm, colors

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
                 'p':   {'Type':'product', 'Formula':'dataPlot * np.conj(dataPlot)',  'Exec': lambda x: x * np.conj(x)},
                 'phase':{'Type':'phase', 'Formula':'np.angle(dataPlot)',  'Exec': lambda x: np.angle(x)},
                 'phaseUW':{'Type':'phase (unwrapped)', 'Formula':'np.unwrap(np.angle(dataPlot), axis = axisUW)',  'Exec': lambda x: np.unwrap(np.angle(x), axis = axisNum)}
                 }

    if returnDict:
        return pTypeDict

    # Set pType in output plotting data Xarray
    if pType is 'phaseUW':
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
def sphFromBLMPlot(BLMXin, res = 50, pType = 'r', plotFlag = False, facetDim = None, backend = 'mpl', fnType = None):
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

    Returns
    -------
    Xarray
        Xarray containing the calcualted surfaces I(theta,phi,...)

    fig
        List of figure handles.

    '''

    #*** Check dataType

    # For ADMs, can use this plotting routine if S=0 only.
    # NOTE - this should work if other singleton dims are present, but will need to explicitly set facetDims for plots, otherwise will throw an error
    if (BLMXin.attrs['dataType'] is 'ADM'):
        if all(BLMXin.S == 0):
            # NOTE: Squeeze here to kill S dim, but this will also drop (l,m) if singleton. Needs a fix.
            # TODO: fix this!
            BLMX = BLMXin.copy().unstack('ADM').rename({'K':'l', 'Q':'m'}).drop('S').stack({'BLM':('l','m')}).squeeze()
        else:
            print('***ADM set with non-zero S, skipping Ylm routine.')
            BLMX = None

    # For normal BLMs, no action required.
    if BLMXin.attrs['dataType'] is 'BLM':
        BLMX = BLMXin

    fig = None
    if BLMX is not None:
        # Calculate YLMs
        if hasattr(BLMX,'normType') and (fnType is None):
            fnType = BLMX.attrs['normType']
            print(f'Using {fnType} betas (from BLMX array).')
        elif fnType is not None:
            print(f'Using {fnType} betas (as passed).')
        else:
            fnType = 'sph'
            print(f'Using default {fnType} betas.')


        YLMX = sphCalc(BLMX.l.max(), res=res, fnType=fnType)
        YLMX = YLMX.rename({'LM':'BLM'})    # Switch naming for multiplication & plotting

        # Calculate MFPADs (theta,phi)
        dataPlot = BLMX*YLMX
        dataPlot = dataPlot.rename({'BLM':'LM'})    # Switch naming back for plotting function
        dataPlot.attrs['normType'] = fnType

        # Pass data to plotting function
        if plotFlag:
            fig = sphSumPlotX(dataPlot, pType = pType, facetDim = facetDim, backend = backend)


    return dataPlot, fig


# Sum and plot spherical functions from an Xarray.
# TODO: This currently assumes matplotlib cm is already loaded.
# TODO: More plot types.
# TODO: More careful testing - not totally sure if summation & axes consistent here.
def sphSumPlotX(dataIn, pType = 'a', facetDim = 'Eke', backend = 'mpl',  convention = 'phys'):
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

    backend : str, optional, default 'mpl'
        Set backend used for plotting.

        - mpl   matplotlib: basic 3D plotting, one figure per surface.
        - pl    plotly: fancier 3D plotting, interactive in Jupyter but may fail at console.
                Subplots for surfaces.
        - hv    holoviews: fancier plotting with additional back-end options.
                Can facet on specific data types.

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


    '''

    # Sum over QNs if necessary
    # if (dataIn.LM.shape[0]) > 1:
    if 'LM' in dataIn.dims:
        dataPlot = dataIn.sum(dim='LM')
    else:
        dataPlot = dataIn

    # (theta,phi) grid from Xarray coords
    theta, phi = np.meshgrid(dataPlot.Theta, dataPlot.Phi)

    # Set data according to type of plot selected
    dataPlot = plotTypeSelector(dataPlot, pType = pType)

    # Switch plot function based on backend
    print('Plotting with {0}'.format(backend))

    fig = []

    # Matplotlib
    if backend is 'mpl':
        # Check dimensionality - loop over facetDim if necessary
        # Return list of fig handles
        if len(dataPlot.dims) > 2:
            print('Data dims: {0}, subplots on {1}'.format(dataPlot.dims, facetDim))

            # Basic loop over facetDim
            # Fails for dims with multi-index, or repeated values for sub-indexes if used as selector.
            # Therefore, use positional numerical index...
            # for nPlot in dataPlot[facetDim]:
            #     fig.append(sphPlotMPL(dataPlot.sel({facetDim:[nPlot]}).squeeze(), theta, phi))  # Will fail if dims>2 passed.

            # Loop over facetDim with full dataArray as selector, allows for MultiIndex cases.
            for nPlot in dataPlot[facetDim]:
                fig.append(sphPlotMPL(dataPlot.sel({facetDim:nPlot}).squeeze(), theta, phi,  convention = convention))  # Will fail if dims>2 passed.

        else:
            # Call matplotlib plotting fn., single surface
            fig.append(sphPlotMPL(dataPlot, theta, phi,  convention = convention))

    # Plotly
    if backend is 'pl':
        fig.append(sphPlotPL(dataPlot, theta, phi, facetDim))

    return fig


#******** Low-level plotting functions

# Set cart coords from spherical polar coords
def sphToCart(R, theta, phi, convention = 'phys'):
    r"""
    Convert spherical polar coords :math:`(R,\theta,\phi)` to Cartesian :math:`(X,Y,Z)`.

    Parameters
    ----------
    R, theta, phi : np.arrays
        Spherical polar coords :math:`(R,\theta,\phi)`.

    convention : str, optional, default = 'phys'
        Specify choice of Spherical Polar coordinate system, 'phys' or 'maths' (see note below).

    Returns
    -------
    X, Y, Z : np.arrays
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

    return X, Y, Z


# Plot using matplotlib - legacy
# TODO - remove this.
# def sphPlot(dataPlot, theta, phi):
#     sphPlotMatplotLib(dataPlot, theta, phi)

# Plot using matplotlib
def sphPlotMPL(dataPlot, theta, phi, convention = 'phys'):
    '''
    Plot spherical polar function (R,theta,phi) to a Cartesian grid, using Matplotlib.

    Parameters
    ----------
    dataPlot : np.array or Xarray
        Values to plot, single surface only, with dims (theta,phi).

    theta, phi : np.arrays
        Angles defining spherical polar grid, 2D arrays.

    Returns
    -------
    fig
        Handle to matplotlib figure.


    '''

    X, Y, Z = sphToCart(dataPlot, theta, phi, convention = convention)

    # Plot in a new figure using Matplotlib
    fig = plt.figure()
    ax = Axes3D(fig)
    # ax.plot_surface(x, y, z, rstride=1, cstride=1, cmap=cm.jet)
    ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap=cm.jet)
    # ax.axis('equal') # Not implemented for 3D axes
    # Rescale axis to equal, hack from https://stackoverflow.com/questions/8130823/set-matplotlib-3d-plot-aspect-ratio
    scaling = np.array([getattr(ax, 'get_{}lim'.format(dim))() for dim in 'xyz'])
    ax.auto_scale_xyz(*[[np.min(scaling), np.max(scaling)]]*3)
    plt.show()

    return fig


# Plot as Plotly subplots, as function of selected dim.
# Currently set for subplotting over facetDim, but assumes other dims (theta,phi)
def sphPlotPL(dataPlot, theta, phi, facetDim = 'Eke', rc = None):
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

    Returns
    -------
    fig
        Handle to figure.

    '''

    # Set up subplots
    nData = dataPlot[facetDim].size

    if rc is None:
        nCols = 3
        rc = [nData/nCols, nCols]

    pType = {'type':'surface'}
    specs = [[pType] * rc[1] for i in range(rc[0])]  # Set specs as 2D list of dicts.

    fig = make_subplots(rows=rc[0], cols=rc[1], specs=specs)
    nPlots = rc[0] * rc[1]

    # Add surfaces
    # From https://plot.ly/python/3d-subplots/

    n=0
    for rInd in range(1,rc[0]+1):
        for cInd in range(1,rc[1]+1):
            if n < nData:
                X,Y,Z = sphToCart(dataPlot.sel({facetDim:dataPlot[facetDim][n]}),theta,phi)  # Bit ugly, probably a better way to select here.
                fig.add_trace(
                    go.Surface(x=X, y=Y, z=Z, colorscale='Viridis', showscale=False),
                    row=rInd, col=cInd) # , title=[facetDim + '=' + str(dataPlot[facetDim][n].data)])  # Needs some work here...

                # fig['layout'].update(scene=dict(aspectmode="data"))  # Try and fix aspect ratio issues - getting stuck on first subplot?  Seems independent of setting here.

                # DOESN'T WORK
                # From https://github.com/plotly/plotly.py/issues/70
                # fig['layout'].update(scene=dict(
                #                     aspectmode='manual',
                #                     aspectratio=go.layout.scene.Aspectratio(
                #                     x=x.ptp(), y=y.ptp(), z=z.pyp())))

                n=n+1

    fig.show()

    return fig
