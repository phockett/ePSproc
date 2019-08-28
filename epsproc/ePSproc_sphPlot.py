# -*- coding: utf-8 -*-
"""
ePSproc spherical polar plotting functions

Collection of functions for plotting I(theta, phi) type data.

14/08/19    v1

"""

# imports
import numpy as np
# For plotting functions
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm, colors

# Sum and plot spherical functions from an Xarray.
# TODO: This currently assumes matplotlib cm is already loaded.
# TODO: More plot types.
# TODO: More careful testing - not totally sure if summation & axes consistent here.
def sphSumPlotX(dataIn, pType = 'a', backend = 'matplotlib'):
    '''
    Plot sum of spherical harmonics from an Xarray.

    Parameters
    ----------
    dataIn : Xarray
        Input set of precalculated Ylms, dims (theta,phi) or (lm,theta,phi).
        In the latter case (lm) dimension is summed.

    pType : char, optional, default 'a'
        Set type of plot.
        'a' = np.abs(Ylm)
        'r' = np.real(Ylm)
        'i' = np.imag(Ylm)
        'p' = Ylm * np.conj(Ylm)

    backend : str, optional, default 'matplotlib'
        Set backend used for plotting.
        - matplotlib: basic 3D plotting, one figure per surface.
        - plotly: fancier 3D plotting, interactive in Jupyter but may fail at console.
                    Subplots for surfaces.
        - holoviews: fancier plotting with additional back-end options.
                    Can facet on specific data types.

    Example
    -------
    >>> YlmX = sphCalc(2, res = 50)
    >>> sphSumPlotX(YlmX)

    Notes
    -----
    Pretty basic functionality here, should add more colour mapping options and multiple plots.

    '''

    # Sum over QNs if necessary
    # if (dataIn.LM.shape[0]) > 1:
    if 'LM' in dataIn.dims:
        dataPlot = dataIn.sum(dim='LM')
    else:
        dataPlot = dataIn

    # (theta,phi) grid from Xarray coords
    theta, phi = np.meshgrid(dataPlot.Theta, dataPlot.Phi)

    # Check dimensionality - loop over dims if necessary
    if len(dataPlot.dims) > 2:
        print('Data dims: {0}'.format(dataPlot.dims))

    else:
        # Call basic plotting fn., single surface
        sphPlot(dataPlot, theta, phi, pType, backend)


#******** Low-level plotting functions

# Set cart coords from spherical polar coords
def sphToCart(R, theta, phi):
    """
    Convert spherical polar coords (R,theta,phi) to Cartesian (X,Y,Z).

    Parameters
    ----------
    R, theta, phi : np.arrays defining a grid of Sph coords.

    Outputs
    -------
    X, Y, Z : np.arrays defining a grid of Cart coords.

    """

    X = R * np.sin(phi) * np.cos(theta)
    Y = R * np.sin(phi) * np.sin(theta)
    Z = R * np.cos(phi)

    return X, Y, Z

def sphPlot(dataPlot, theta, phi, pType = 'a', backend = 'matplotlib'):
    '''
    Plot spherical polar function (R,theta,phi) to a Cartesian grid, using Matplotlib.

    Parameters
    ----------
    dataPlot : np.array
        Values to plot, dims (theta,phi).

    theta, phi: np.arrays
        Angles defining spherical polar grid, 2D arrays.

    pType : char, optional, default 'a'
        Set type of plot.
        'a' = np.abs(dataPlot)
        'r' = np.real(dataPlot)
        'i' = np.imag(dataPlot)
        'p' = dataPlot * np.conj(dataPlot)

    Example
    -------
    >>> YlmX = sphCalc(2, res = 50)
    >>> sphSumPlotX(YlmX)

    Notes
    -----
    Pretty basic functionality here, should add more colour mapping options and multiple plots, and alternative back-ends (Matplotlib only currently).
    Some functionality may end up being pushed to plotting library, e.g. Holoviews, which already supports faceting etc. etc.

    '''

    # Convert to cart
    if pType == 'a':
        R = np.abs(dataPlot)
    elif pType == 'r':
        R = np.real(dataPlot)
    elif pType == 'i':
        R = np.imag(dataPlot)
    elif pType == 'p':
        # TODO: check conj defns with Xarray
        R = np.abs(dataPlot * np.conj(dataPlot))

    X, Y, Z = sphToCart(R,theta,phi)

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
