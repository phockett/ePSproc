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
    elif BLMXin.attrs['dataType'] is 'BLM':
        BLMX = BLMXin

    # Default case - set as BLM case but include dataType warning.
    else:
        print(f"*** sphPlot dataType = {BLMXin.attrs['dataType']} not recognised, trying anyway.")
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
        dataPlot.attrs = BLMX.attrs # Ensure attrs propagated
        dataPlot.attrs['normType'] = fnType

        # Pass data to plotting function
        if plotFlag:
            fig = sphSumPlotX(dataPlot, pType = pType, facetDim = facetDim, backend = backend)


    return dataPlot, fig


# Sum and plot spherical functions from an Xarray.
# TODO: This currently assumes matplotlib cm is already loaded.
# TODO: More plot types.
# TODO: More careful testing - not totally sure if summation & axes consistent here.
def sphSumPlotX(dataIn, pType = 'a', facetDim = 'Eke', backend = 'mpl',  convention = 'phys', titleString = None):
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

    convention : str, optional, default = 'phys'
        Spherical polar coord convention, see :py:func:`epsproc.sphToCart`

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

    # Set main plot title if not passed.
    if titleString is None:
        if hasattr(dataIn,'jobLabel'):
            titleString = dataIn.attrs['jobLabel']
        elif hasattr(dataIn,'file'):
            titleString = dataIn.attrs['file']
        else:
            titleString = ""

    print(f"Sph plots: {titleString}")  # Throw out titlestring here for ease.


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
                # fig.append(sphPlotMPL(dataPlot.sel({facetDim:nPlot}).squeeze(), theta, phi,  convention = convention, tString = f"{facetDim}: {nPlot.item()}"))
                # {dataPlot[facetDim][nPlot].item()}"))  # Will fail if dims>2 passed.
                fig.append(sphPlotMPL(dataPlot.sel({facetDim:nPlot}).squeeze(), theta, phi,  convention = convention, tString = titleString + f"\n{facetDim}: {nPlot.item()}"))

        else:
            # Call matplotlib plotting fn., single surface
            fig.append(sphPlotMPL(dataPlot, theta, phi,  convention = convention))

    # Plotly - note that faceting is handled directly by Plotly in this case.
    if backend is 'pl':
        fig.append(sphPlotPL(dataPlot, theta, phi, facetDim, convention = convention))

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
def sphPlotMPL(dataPlot, theta, phi, convention = 'phys', tString = None):
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
def sphPlotPL(dataPlot, theta, phi, facetDim = 'Eke', rc = None, norm = 'global', convention = 'phys'):
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

    rc : array, optional, default = None
        If set, use to define layout grid [rows, columns].
        If not set, this will be set for nCols = 5

    norm : str, optional, default = 'global'
        Set for plot normalisation.
        - 'global' : use same (x,y,z) limits for all plots.
        - 'local' : auto set (x,y,z) limits for each plot.

    convention : str, optional, default = 'phys'
        Spherical polar coord convention, see :py:func:`epsproc.sphToCart`


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

    '''

    # Set up subplots
    nData = dataPlot[facetDim].size

    if rc is None:
        nCols = 6
        if nData < nCols:   # Set for single row layout
            nCols = nData
        rc = [nData/nCols, nCols]

    # rc = np.round(rc).astype(np.int)  # Fix dtype - Plotly throws type error for np types however.
    rc = [int(np.ceil(rc[0])), int(np.ceil(rc[1]))]

    pType = {'type':'surface'}
    specs = [[pType] * rc[1] for i in range(rc[0])]  # Set specs as 2D list of dicts.

    # Check/set some global props - assumes Xarray datatype at the moment (also assumed later!)
    if isinstance(dataPlot, xr.core.dataarray.DataArray):
        # Check max coord limit, to use for all subplot ranges
        # Should rewrite sphToCart for Xarray output?
        X,Y,Z =sphToCart(dataPlot,dataPlot.Theta,dataPlot.Phi)

        Cmax = np.max([X.max(), Y.max(), Z.max()])  # Check max & min coord values - may be unequal in some cases.
        Cmin = np.min([X.min(), Y.min(), Z.min()])
        Rmax = np.max([np.abs(Cmin), Cmax])
        padding = 0.1 * Rmax  # Padding %age
        aRanges = dict(range=[-(Rmax + padding), Rmax+padding])

        # Set subplot titles, see https://stackoverflow.com/questions/53991810/how-to-add-subplot-title-to-3d-subplot-in-plotly
        titles = [f"{facetDim}: {item.item()}" for item in dataPlot[facetDim]]

    # Set up subplots
    fig = make_subplots(rows=rc[0], cols=rc[1], specs=specs, subplot_titles=titles)
    nPlots = rc[0] * rc[1]

    # Add surfaces
    # From https://plot.ly/python/3d-subplots/

    n=0
    for rInd in range(1,rc[0]+1):
        for cInd in range(1,rc[1]+1):
            if n < nData:

                # Set data (NOTE - now repeats above in Xarray case)
                X,Y,Z = sphToCart(dataPlot.sel({facetDim:dataPlot[facetDim][n]}),theta,phi)  # Bit ugly, probably a better way to select here.

                # Set plot object
                showscale = False
                # if rInd == (rc[0]):
                #     showscale = True  # Set scale bar for row?  Only for global scaling? With this code get multiple cbars?

                fig.add_trace(
                    go.Surface(x=X, y=Y, z=Z, colorscale='Viridis', showscale=showscale),
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
                if norm is 'global':
                    # Try looping... OK with dict unpacking... huzzah!
                    # NOTE Scene indexing starts at 1, so do this after n increments
                    options = dict(xaxis = aRanges, yaxis = aRanges, zaxis = aRanges, aspectmode='cube')

                else:
                    options = dict(aspectmode='cube')

                fig.update_layout(**{sceneN:options})  # No effect of aspect here? auto/cube/data/manual

                # if rInd == (rc[0]):
                #     fig.update_layout(**{sceneN:dict(showscale=True)})  # Set scale bar for row?


    fig.show()

    return fig
