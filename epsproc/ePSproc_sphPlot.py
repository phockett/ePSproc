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
def sphSumPlotX(dataIn, pType = 'a'):
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

    # Call basic plotting fn.
    sphPlot(dataPlot, theta, phi, pType)

    # THIS CODE NOW MOVED TO sphPlot()
#    # Convert to cart
#    if pType == 'a':
#        R = np.abs(dataPlot)
#    elif pType == 'r':
#        R = np.real(dataPlot)
#    elif pType == 'i':
#        R = np.imag(dataPlot)
#    elif pType == 'p':
#        # TODO: check conj defns with Xarray
#        R = np.abs(dataPlot * np.conj(dataPlot))
#
#    X = R * np.sin(phi) * np.cos(theta)
#    Y = R * np.sin(phi) * np.sin(theta)
#    Z = R * np.cos(phi)
#
#    # Plot in a new figure using Matplotlib
#    fig = plt.figure()
#    ax = Axes3D(fig)
#    # ax.plot_surface(x, y, z, rstride=1, cstride=1, cmap=cm.jet)
#    ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap=cm.jet)
#    plt.show()

# Low-level sph plotter
def sphPlot(dataPlot, theta, phi, pType = 'a'):
    '''
    Plot spherical polar function (R,theta,phi) to a Cartesian grid.

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

    X = R * np.sin(phi) * np.cos(theta)
    Y = R * np.sin(phi) * np.sin(theta)
    Z = R * np.cos(phi)

    # Plot in a new figure using Matplotlib
    fig = plt.figure()
    ax = Axes3D(fig)
    # ax.plot_surface(x, y, z, rstride=1, cstride=1, cmap=cm.jet)
    ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap=cm.jet)
    plt.show()
