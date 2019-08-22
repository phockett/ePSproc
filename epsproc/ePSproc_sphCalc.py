# -*- coding: utf-8 -*-
"""
ePSproc spherical function calculations.

Collection of functions for calculating Ylm, wignerD etc.

14/08/19    v1

"""

# Imports
import numpy as np
import pandas as pd
import xarray as xr
from scipy.special import sph_harm

# Calculate a set of sph function
def sphCalc(Lmax, Lmin = 0, res = None, angs = None, XFlag = True):
    '''
    Calculate set of spherical harmonics Ylm(theta,phi) on a grid.

    Parameters
    ----------
    Lmax : int
        Maximum L for the set. Ylm calculated for Lmin:Lmax, all m.
    Lmin : int, optional, default 0
        Min L for the set. Ylm calculated for Lmin:Lmax, all m.
    res : int, optional, default None
        (Theta, Phi) grid resolution, outputs will be of dim [res,res].
    angs : list of 2D np.arrays, [thetea, phi], optional, default None
        If passed, use these grids for calculation
    XFlag : bool, optional, default True
        Flag for output. If true, output is Xarray. If false, np.arrays

    Note that either res OR angs needs to be passed.

    Outputs
    -------
    - if XFlag -
    YlmX
        3D Xarray, dims (lm,theta,phi)
    - else -
    Ylm, lm
        3D np.array of values, dims (lm,theta,phi), plus list of lm pairs

    Methods
    -------
    Currently set for scipy.special.sph_harm as calculation routine.

    Example
    -------
    >>> YlmX = sphCalc(2, res = 50)

    '''

    # Set coords based on inputs
    # TODO: better code here (try/fail?)
    if angs is None and res:
        theta, phi = np.meshgrid(np.linspace(0,2*np.pi,res),np.linspace(0,np.pi,res))
    elif res is None and angs:
        theta = angs[0]
        phi = angs[1]
    else:
        print('Need to pass either res or angs.')
        return False

    # Loop over lm and calculate
    lm = []
    Ylm = []
    for l in np.arange(Lmin,Lmax+1):
        for m in np.arange(-l,l+1):
            lm.append([l, m])
            Ylm.append(sph_harm(m,l,theta,phi))

    # Return as Xarray or np arrays.
    if XFlag:
        # Set indexes
        QNs = pd.MultiIndex.from_arrays(np.asarray(lm).T, names = ['l','m'])
        YlmX = xr.DataArray(np.asarray(Ylm), coords=[('QN',QNs), ('Theta',theta[0,:]), ('Phi',phi[:,0])])
        return YlmX
    else:
        return np.asarray(Ylm), np.asarray(lm)
