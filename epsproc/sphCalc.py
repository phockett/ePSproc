# -*- coding: utf-8 -*-
"""
ePSproc spherical function calculations.

Collection of functions for calculating Spherical Tensors: Ylm, wignerD etc.

For spherical harmonics, currently using scipy.special.sph_harm

For other functions, using Moble's spherical_functions package
https://github.com/moble/spherical_functions

See tests/Spherical function testing Aug 2019.ipynb

02/12/19        Added basic TKQ multipole frame rotation routine.
27/08/19        Added wDcalc for Wigner D functions.
14/08/19    v1  Implmented sphCalc

"""

# Imports
import numpy as np
import pandas as pd
import xarray as xr
from scipy.special import sph_harm
import spherical_functions as sf
import quaternion
from sympy.physics.quantum.spin import Rotation  # For basic frame rotation code, should update to use sf

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
        YlmX = xr.DataArray(np.asarray(Ylm), coords=[('LM',QNs), ('Theta',theta[0,:]), ('Phi',phi[:,0])])
        return YlmX
    else:
        return np.asarray(Ylm), np.asarray(lm)


# Calculate wignerD functions
#   Adapted directly from Matlab code,
#   via Jupyter test Notebook "Spherical function testing Aug 2019.ipynb"
def wDcalc(Lrange = [0, 1], Nangs = None, eAngs = None, R = None, XFlag = True):
    '''
    Calculate set of Wigner D functions D(l,m,mp,R) on a grid.

    Parameters
    ----------
    Lrange : list, optional, default [0, 1]
        Range of L to calculate parameters for.
        If len(Lrange) == 2 assumed to be of form [Lmin, Lmax], otherwise list is used directly.
        For a given l, all (m, mp) combinations are calculated.

    Options for setting angles (use one only):
    Nangs : int, optional, default None
        If passed, use this to define Euler angles sampled.
        Ranges will be set as (theta, phi, chi) = (0:pi, 0:pi/2, 0:pi) in Nangs steps.
    eAngs : np.array, optional, default None
        If passed, use this to define Euler angles sampled.
        Array of angles, [theta,phi,chi], in radians
    R : np.array, optional, default None
        If passed, use this to define Euler angles sampled.
        Array of quaternions, as given by quaternion.from_euler_angles(eAngs).


    XFlag : bool, optional, default True
        Flag for output. If true, output is Xarray. If false, np.arrays


    Outputs
    -------
    - if XFlag -
    wDX
        Xarray, dims (lmmp,Euler)
    - else -
    wD, R, lmmp
        np.arrays of values, dims (lmmp,Euler), plus list of angles and lmmp sets.

    Methods
    -------
    Uses Moble's spherical_functions package for wigner D function.
    https://github.com/moble/spherical_functions

    Moble's quaternion package for angles and conversions.
    https://github.com/moble/quaternion

    Examples
    --------
    >>> wDX1 = wDcalc(eAngs = np.array([0,0,0]))

    >>> wDX2 = wDcalc(Nangs = 10)

    '''
    # Set QNs for calculation, (l,m,mp)
    if len(Lrange) == 2:
        Ls = np.arange(Lrange[0], Lrange[1]+1)
    else:
        Ls = Lrange

    QNs = []

    for l in Ls:
        for m in np.arange(-l, l+1):
            for mp in np.arange(-l, l+1):
                QNs.append([l, m, mp])

    QNs = np.array(QNs)

    # Set angles - either input as a range, a set or as quaternions
    if Nangs is not None:
        # Set a range of Eugler angles for testing
        pRot = np.linspace(0,np.pi,Nangs)
        tRot = np.linspace(0,np.pi/2,Nangs)
        cRot = np.linspace(0,np.pi,Nangs)
        eAngs = np.array([pRot, tRot, cRot,]).T

    if eAngs is not None:
        if eAngs.shape[-1] != 3:    # Check dims, should be (N X 3) for quaternion... but transpose for pd.MultiIndex
            eAngs = eAngs.T

    if R is None:
        # Convert to quaternions
        R =  quaternion.from_euler_angles(eAngs)


    # Calculate WignerDs
    # sf.Wigner_D_element is vectorised for QN OR angles
    # Here loop over QNs for a set of angles R
    wD = []
    lmmp = []
    for n in np.arange(0, QNs.shape[0]):
        lmmp.append(QNs[n,:])
        wD.append(sf.Wigner_D_element(R, QNs[n,0], QNs[n,1], QNs[n,2]))

    # Return values as Xarray or np.arrays
    if XFlag:
        # Put into Xarray
        #TODO: this will currently fail for a single set of QNs.
        QNs = pd.MultiIndex.from_arrays(np.asarray(lmmp).T, names = ['lp','mu','mu0'])
        if (eAngs is not None) and (eAngs.size == 3):  # Ugh, special case for only one set of angles.
            Euler = pd.MultiIndex.from_arrays([[eAngs[0]],[eAngs[1]],[eAngs[2]]], names = ['P','T','C'])
            wDX = xr.DataArray(np.asarray(wD), coords=[('QN',QNs)])
            wDX = wDX.expand_dims({'Euler':Euler})
        else:
            Euler = pd.MultiIndex.from_arrays(eAngs.T, names = ['P','T','C'])
            wDX = xr.DataArray(np.asarray(wD), coords=[('QN',QNs), ('Euler',Euler)])

        return wDX

    else:
        return wD, R, np.asarray(lmmp).T



#*** Basic frame rotation code, see https://github.com/phockett/Quantum-Metrology-with-Photoelectrons/blob/master/Alignment/Alignment-1.ipynb
# Define frame rotation of state multipoles.
# Eqn. 4.41 in Blum (p127)
# Currently a bit ugly!
# Also set for numerical output only, although uses Sympy functions which can be used symbolically.
# Pass TKQ np.array [K,Q,TKQ], eAngs list of Euler angles (theta,phi,chi) to define rotation.
def TKQarrayRot(TKQ,eAngs):
    r"""
    Frame rotation for multipoles $T_{K,Q}$.

    Basic frame rotation code, see https://github.com/phockett/Quantum-Metrology-with-Photoelectrons/blob/master/Alignment/Alignment-1.ipynb for examples.

    Parameters
    ----------
    TKQ : np.array
        Values defining the initial distribution, [K,Q,TKQ]

    eAngs : list or np.array
        List of Euler angles (theta,phi,chi) defining rotated frame.

    Returns
    -------

    TKQRot : np.array
        Multipoles $T'_{K,Q}$ in rotated frame, as an np.array [K,Q,TKQ].

    TODO: redo with Moble's functions, and Xarray input & output.

    Formalism
    ----------

    For the state multipoles, frame rotations are fairly straightforward
    (Eqn. 4.41 in Blum):

    .. math::
        \begin{equation}
        \left\langle T(J',J)_{KQ}^{\dagger}\right\rangle =\sum_{q}\left\langle T(J',J)_{Kq}^{\dagger}\right\rangle D(\Omega)_{qQ}^{K*}
        \end{equation}

    Where $D(\Omega)_{qQ}^{K*}$ is a Wigner rotation operator, for a
    rotation defined by a set of Euler angles $\Omega=\{\theta,\phi,\chi\}$.
    Hence the multipoles transform, as expected, as irreducible tensors,
    i.e. components $q$ are mixed by rotation, but terms of different
    rank $K$ are not.

    """

    TKQRot = []
    thres = 1E-5
    Kmax = 6

    # Easy way - loop over possible output values & sum based on input TKQ. Can probably do this in a smarter way.
    for K in range(0,Kmax+1):
        for q in range(-K,K+1):

            # Set summation variable and add relevant terms from summation
            TKQRotSum = 0.0
            for row in range(TKQ.shape[0]):
                Kin = TKQ[row][0]
                Qin = TKQ[row][1]

                if Kin == K:
                    Dval = Rotation.D(K,Qin,q,eAngs[0],eAngs[1],eAngs[2])
                    TKQRotSum += conjugate(Dval.doit())*TKQ[row][2]
                else:
                    pass

            if np.abs(N(TKQRotSum)) > thres:
                TKQRot.append([K,q,N(TKQRotSum)])  # Use N() here to ensure Sympy numerical output only

    return np.array(TKQRot)
