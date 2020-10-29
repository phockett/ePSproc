# -*- coding: utf-8 -*-
r"""
ePSproc MFPAD functions
-----------------------

29/10/20    v2  Updated with Wigner delay routine.
                Tidied up To Do list items.
                Tidied up docstring.

05/08/19    v1  Initial python version.
                Based on original Matlab code ePS_MFPAD.m

Structure
---------
Calculate MFPAD on a grid from input ePS matrix elements.
Use fast functions, pre-calculate if possible.
Data in Xarray, use selection functions and multiplications based on relevant quantum numbers, other axes summed over.

Choices for functions...
    * `Moble's spherical functions (quaternion based) <https://github.com/moble/spherical_functions>`_

      Provides fast wignerD, 3j and Ylm functions, uses Numba.

      Install with:

      >>> conda install -c conda-forge spherical_functions

    * `Scipy special.sph_harm <https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.sph_harm.html#scipy.special.sph_harm>`_



To Do
-----
* More efficient computation?  Use Xarray group by?

Formalism
---------

MF-PADs numerical expansion
============================

This is implemented numerically in :py:func:`epsproc.MFPAD.mfpad()`. For direct computation of $\beta_{LM}$ parameters (rather than full observable expansion) see :py:func:`epsproc.MFBLM` (original looped version) or :py:func:`epsproc.geomFunc.mfblmGeom` (geometric tensor version).

From `ePSproc: Post-processing suite for ePolyScat electron-molecule scattering calculations <https://www.authorea.com/users/71114/articles/122402/_show_article>`_.

.. math::
    I_{\mu_{0}}(\theta_{\hat{k}},\phi_{\hat{k}},\theta_{\hat{n}},\phi_{\hat{n}})=\frac{4\pi^{2}E}{cg_{p_{i}}}\sum_{\mu_{i},\mu_{f}}|T_{\mu_{0}}^{p_{i}\mu_{i},p_{f}\mu_{f}}(\theta_{\hat{k}},\phi_{\hat{k}},\theta_{\hat{n}},\phi_{\hat{n}})|^{2}\label{eq:MFPAD}

    T_{\mu_{0}}^{p_{i}\mu_{i},p_{f}\mu_{f}}(\theta_{\hat{k}},\phi_{\hat{k}},\theta_{\hat{n}},\phi_{\hat{n}})=\sum_{l,m,\mu}I_{l,m,\mu}^{p_{i}\mu_{i},p_{f}\mu_{f}}(E)Y_{lm}^{*}(\theta_{\hat{k}},\phi_{\hat{k}})D_{\mu,\mu_{0}}^{1}(R_{\hat{n}})\label{eq:TMF}

    I_{l,m,\mu}^{p_{i}\mu_{i},p_{f}\mu_{f}}(E)=\langle\Psi_{i}^{p_{i},\mu_{i}}|\hat{d_{\mu}}|\Psi_{f}^{p_{f},\mu_{f}}\varphi_{klm}^{(-)}\rangle\label{eq:I}


In this formalism:

* :math:`I_{l,m,\mu}^{p_{i}\mu_{i},p_{f}\mu_{f}}(E)` is the radial part of the dipole matrix element, determined from the initial and final state electronic wavefunctions :math:`\Psi_{i}^{p_{i},\mu_{i}}` and :math:`\Psi_{f}^{p_{f},\mu_{f}}`, photoelectron wavefunction :math:`\varphi_{klm}^{(-)}` and dipole operator :math:`\hat{d_{\mu}}`. Here the wavefunctions are indexed by irreducible representation (i.e. symmetry) by the labels :math:`p_{i}` and :math:`p_{f}`, with components :math:`\mu_{i}` and :math:`\mu_{f}` respectively; :math:`l,m` are angular momentum components, :math:`\mu` is the projection of the polarization into the MF (from a value :math:`\mu_{0}` in the LF). Each energy and irreducible representation corresponds to a calculation in ePolyScat.
* :math:`T_{\mu_{0}}^{p_{i}\mu_{i},p_{f}\mu_{f}}(\theta_{\hat{k}},\phi_{\hat{k}},\theta_{\hat{n}},\phi_{\hat{n}})` is the full matrix element (expanded in polar coordinates) in the MF, where :math:`\hat{k}` denotes the direction of the photoelectron :math:`\mathbf{k}`-vector, and :math:`\hat{n}` the direction of the polarization vector :math:`\mathbf{n}` of the ionizing light. Note that the summation over components :math:`\{l,m,\mu\}` is coherent, and hence phase sensitive.
* :math:`Y_{lm}^{*}(\theta_{\hat{k}},\phi_{\hat{k}})` is a spherical harmonic.
* :math:`D_{\mu,\mu_{0}}^{1}(R_{\hat{n}})` is a Wigner rotation matrix element, with a set of Euler angles :math:`R_{\hat{n}}=(\phi_{\hat{n}},\theta_{\hat{n}},\chi_{\hat{n}})`, which rotates/projects the polarization into the MF.
* :math:`I_{\mu_{0}}(\theta_{\hat{k}},\phi_{\hat{k}},\theta_{\hat{n}},\phi_{\hat{n}})` is the final (observable) MFPAD, for a polarization :math:`\mu_{0}` and summed over all symmetry components of the initial and final states, :math:`\mu_{i}` and :math:`\mu_{f}`. Note that this sum can be expressed as an incoherent summation, since these components are (by definition) orthogonal.
* :math:`g_{p_{i}}` is the degeneracy of the state :math:`p_{i}`.

For more details see: Toffoli, D., Lucchese, R. R., Lebech, M., Houver, J. C., & Dowek, D. (2007). Molecular frame and recoil frame photoelectron angular distributions from dissociative photoionization of NO2. The Journal of Chemical Physics, 126(5), 054307. https://doi.org/10.1063/1.2432124


(MF) Wigner Delays
==================

This is implemented numerically in :py:func:`epsproc.MFPAD.mfWignerDelay()`.

The partial wave, or channel resolved, Wigner delay can be given as the energy-derivative of the continuum wavefunction phase, $\arg(\psi_{lm}^{*})$

.. math::
    \begin{equation}
    \tau_{w}(k,\theta,\phi)=\hbar\frac{d\arg(\psi_{lm}^{*})}{d\epsilon}
    \end{equation}


And the total, or group, delay from the (coherent) sum over these channels.

.. math::
    \begin{equation}
    \tau_{w}^{g}(k,\theta,\phi)=\hbar\frac{d\arg(\sum_{l,m}\psi_{lm}^{*})}{d\epsilon}\label{eq:tauW_mol_sum}
    \end{equation}


Here the full wavefunction $\sum_{l,m}\psi_{lm}^{*}$ is essentially equivalent to the phase of the (conjugate) $T_{\mu_{0}}^{p_{i}\mu_{i},p_{f}\mu_{f}*}(\theta_{\hat{k}},\phi_{\hat{k}},\theta_{\hat{n}},\phi_{\hat{n}})$ functions given above (and by :py:func:`epsproc.MFPAD.mfpad()`).

For more details see, for example:

- Hockett, Paul, Eugene Frumker, David M Villeneuve, and Paul B Corkum. “Time Delay in Molecular Photoionization.” Journal of Physics B: Atomic, Molecular and Optical Physics 49, no. 9 (May 2016): 095602. https://doi.org/10.1088/0953-4075/49/9/095602.
- Carvalho, C.A.A. de, and H.M. Nussenzveig. “Time Delay.” Physics Reports 364, no. 2 (June 2002): 83–174. https://doi.org/10.1016/S0370-1573(01)00092-8.
- Wigner, Eugene. “Lower Limit for the Energy Derivative of the Scattering Phase Shift.” Physical Review 98, no. 1 (April 1955): 145–47. https://doi.org/10.1103/PhysRev.98.145.

"""


# Imports
import numpy as np
import pandas as pd
import xarray as xr
# Special functions
# from scipy.special import sph_harm
import spherical_functions as sf
import quaternion

import scipy.constants

# Package fns.
# TODO: tidy this up!
from epsproc.util import matEleSelector
from epsproc.sphCalc import sphCalc, setPolGeoms
from epsproc.sphPlot import plotTypeSelector

def mfpad(dataIn, thres = 1e-2, inds = {'Type':'L','it':1}, res = 50, R = None, p = 0):
    """

    Parameters
    ----------
    dataIn : Xarray
        Contains set(s) of matrix elements to use, as output by epsproc.readMatEle().

    thres : float, optional, default 1e-2
        Threshold value for matrix elements to use in calculation.

    ind : dictionary, optional.
        Used for sub-selection of matrix elements from Xarrays.
        Default set for length gauage, single it component only, inds = {'Type':'L','it':'1'}.

    res : int, optional, default 50
        Resolution for output (theta,phi) grids.

    R : list or Xarray of Euler angles or quaternions, optional.
        Define LF > MF polarization geometry/rotations.
        For default case (R = None), calls :py:func:`epsproc.setPolGeoms()`, where 3 geometries are calculated, corresponding to z-pol, x-pol and y-pol cases.
        Defined by Euler angles (p,t,c) = [0 0 0] for z-pol, [0 pi/2 0] for x-pol, [pi/2 pi/2 0] for y-pol.

    p : int, optional.
        Defines LF polarization state, p = -1...1, default p = 0 (linearly pol light along z-axis).
        TODO: add summation over p for multiple pol states in LF.

    Returns
    -------
    Ta
        Xarray (theta, phi, E, Sym) of MFPADs, summed over (l,m)

    Tlm
        Xarray (theta, phi, E, Sym, lm) of MFPAD components, expanded over all (l,m)

    """

    # Define reduced data from selection over all data
    daRed = matEleSelector(dataIn, thres = thres, inds = inds)

    # Generate spherical harmonics
    Lmax = daRed.l.max()
    YlmX = sphCalc(Lmax, res = res)

    # Reindex to match data (should happen automagically, but not always!)
    # YlmXre = YlmX.reindex_like(daRed)

    # Set rotation angles for LF > MF
    if R is None:
        # Set (x,y,z) projection terms only
        # Nangs = 10
        # pRot = np.linspace(0,180,Nangs)
        # tRot = np.linspace(0,90,Nangs)
        # cRot = np.linspace(0,180,Nangs)
        # eAngs = np.array([pRot, tRot, cRot,])*np.pi/180
        # Convert to quaternions
        # R =  quaternion.from_euler_angles(pRot*np.pi/180, tRot*np.pi/180, cRot*np.pi/180)

        # # Eugler angles for rotation of LF->MF, set as [0 0 0] for z-pol, [0 pi/2 0] for x-pol, [pi/2 pi/2 0] for y-pol
        # pRot = [0, 0, np.pi/2]
        # tRot = [0, np.pi/2, np.pi/2]
        # cRot = [0, 0, 0]
        # eAngs = np.array([pRot, tRot, cRot])   # List form to use later
        # Euler = pd.MultiIndex.from_arrays(eAngs, names = ['P','T','C'])
        #
        # # Convert to quaternions
        # R =  quaternion.from_euler_angles(pRot, tRot, cRot)

        # 02/10/20 - updated to use sphCalc.setPolGeoms()
        R = setPolGeoms()



    #**************** Calculate MFPADs

    Tlm = []
    Ta = []

    # Loop over pol geoms R
    for n, Rcalc in enumerate(R):
        T = []
        # Loop over mu terms and multiply
        for mu in np.arange(-1,2):

            # Set by element replacement (preserves whole structure)
            # daTemp = daRed.copy()   # Set explicit copy for rotation.
            # daTemp.loc[{'mu':mu}].values = daTemp.loc[{'mu':mu}].values * sf.Wigner_D_element(Rcalc, 1, mu, 0).conj()

            # Issues with reindexing to extra coords at the moment, so reindex and multiply for specific mu only
            # daTemp = daTemp.sel({'mu':mu})
            # YlmXre = YlmX.reindex_like(daTemp)
            # T.append(YlmXre.conj() * daTemp)  # Output full (l,m,mu) expansion

            # Set by looping and selection
            daTemp = daRed.sel({'mu':mu}) * sf.Wigner_D_element(Rcalc.item(), 1, mu, 0).conj()
            YlmXre = YlmX.reindex_like(daTemp)
            T.append(YlmXre.conj() * daTemp)  # Output full (l,m,mu) expansion

        # Concat & sum over symmetries
        Ts = xr.combine_nested([T[0], T[1], T[2]], concat_dim=['LM'])

        # Add dims - currently set for Euler angles only.
        # Can't seem to add mutiindex as a single element, so set dummy coord here and replace below.
        Ts = Ts.expand_dims({'Euler':[n]})  # Set as index
        # Ts = Ts.expand_dims({'p':[eAngs[0,n]], 't':[eAngs[1,n]], 'c':[eAngs[2,n]]})

        Tlm.append(Ts)
        Ta.append(Ts.sum(dim = 'LM'))

    TlmX = xr.combine_nested(Tlm, concat_dim=['Euler'])
    TaX = xr.combine_nested(Ta, concat_dim=['Euler'])

    # Assign Euler angles to dummy dim
    # TlmX = TlmX.assign_coords(Euler = Euler)
    # TaX = TaX.assign_coords(Euler = Euler)
    TlmX = TlmX.assign_coords(Euler = R['Euler'])  # Coords from Euler Xarray.
    TaX = TaX.assign_coords(Euler = R['Euler'])

    # Set generic data type as TX
    TlmX.attrs['dataType'] = 'TX'
    TaX.attrs['dataType'] = 'TX'

    return TaX, TlmX  # , Ta, Tlm  # For debug also return lists



def mfWignerDelay(dataIn, pType = 'phaseUW'):
    """
    Compute MF Wigner Delay as phase derivative of continuum wavefunction.

    Parameters
    ----------
    dataIn : Xarray
        Contains set(s) of wavefunctions to use, as output by :py:func:`epsproc.MFPAD.mfpad()`.

    pType : str, optional, default = 'phaseUW'
        Used to set data conversion options, as implemented in :py:func:`epsproc.plotTypeSelector()`
        - 'phase' use np.angle()
        - 'phaseUW' unwapped with np.unwrap(np.angle())

    """

    # Set data
    Tw = dataIn.copy()

    Tw['EJ'] = Tw['Eke']*scipy.constants.e  # Convert to J units

    unitConv = scipy.constants.hbar/1e-18  # Convert to attoseconds
    # TODO: move this to util functions.

    # Use existing functionality to set phases - easier if unwrap required (uses lambdas)
    Tw = plotTypeSelector(Tw.conj(), pType = pType).differentiate('EJ')* unitConv

    # Propagate attrs
    Tw.attrs = dataIn.attrs
    Tw.attrs['dataType'] = 'Wigner'
    Tw.attrs['units'] = 'attosecond'

    return Tw 
