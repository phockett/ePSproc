# -*- coding: utf-8 -*-
"""
ePSproc MFPAD functions.


05/08/19    v1  Initial python version.
                Based on original Matlab code ePS_MFPAD.m

Structure
---------
Calculate MFPAD on a grid from input ePS matrix elements.
Use fast functions, pre-calculate if possible.
Data in Xarray, use selection functions and multiplications based on relevant quantum numbers, other axes summed over.

Choices for functions...
- Moble's spherical functions (quaternion based)
     https://github.com/moble/spherical_functions
     conda install -c conda-forge spherical_functions

     wignerD and Ylm functions, quaternion-based

- Scipy
    Ylm: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.sph_harm.html#scipy.special.sph_harm

TODO
----
- Benchmark on NO2 reference results.
- Test over multiple E.
- Test over multiple R.

Formalism
---------
(From `ePSproc: Post-processing suite for ePolyScat electron-molecule scattering calculations <https://www.authorea.com/users/71114/articles/122402/_show_article>`).

.. math::

\begin{equation}
I_{\mu_{0}}(\theta_{\hat{k}},\phi_{\hat{k}},\theta_{\hat{n}},\phi_{\hat{n}})=\frac{4\pi^{2}E}{cg_{p_{i}}}\sum_{\mu_{i},\mu_{f}}|T_{\mu_{0}}^{p_{i}\mu_{i},p_{f}\mu_{f}}(\theta_{\hat{k}},\phi_{\hat{k}},\theta_{\hat{n}},\phi_{\hat{n}})|^{2}\label{eq:MFPAD}
\end{equation}

\begin{equation}
T_{\mu_{0}}^{p_{i}\mu_{i},p_{f}\mu_{f}}(\theta_{\hat{k}},\phi_{\hat{k}},\theta_{\hat{n}},\phi_{\hat{n}})=\sum_{l,m,\mu}I_{l,m,\mu}^{p_{i}\mu_{i},p_{f}\mu_{f}}(E)Y_{lm}^{*}(\theta_{\hat{k}},\phi_{\hat{k}})D_{\mu,\mu_{0}}^{1}(R_{\hat{n}})\label{eq:TMF}
\end{equation}

\begin{equation}
I_{l,m,\mu}^{p_{i}\mu_{i},p_{f}\mu_{f}}(E)=\langle\Psi_{i}^{p_{i},\mu_{i}}|\hat{d_{\mu}}|\Psi_{f}^{p_{f},\mu_{f}}\varphi_{klm}^{(-)}\rangle\label{eq:I}
\end{equation}

In this formalism:
\begin{itemize}
\item $I_{l,m,\mu}^{p_{i}\mu_{i},p_{f}\mu_{f}}(E)$ is the radial part of
the dipole matrix element, determined from the initial and final state
electronic wavefunctions $\Psi_{i}^{p_{i},\mu_{i}}$and $\Psi_{f}^{p_{f},\mu_{f}}$,
photoelectron wavefunction $\varphi_{klm}^{(-)}$ and dipole operator
$\hat{d_{\mu}}$. Here the wavefunctions are indexed by irreducible
representation (i.e. symmetry) by the labels $p_{i}$ and $p_{f}$,
with components $\mu_{i}$ and $\mu_{f}$ respectively; $l,m$ are
angular momentum components, $\mu$ is the projection of the polarization
into the MF (from a value $\mu_{0}$ in the LF). Each energy and irreducible
representation corresponds to a calculation in ePolyScat.
\item $T_{\mu_{0}}^{p_{i}\mu_{i},p_{f}\mu_{f}}(\theta_{\hat{k}},\phi_{\hat{k}},\theta_{\hat{n}},\phi_{\hat{n}})$
is the full matrix element (expanded in polar coordinates) in the
MF, where $\hat{k}$ denotes the direction of the photoelectron $\mathbf{k}$-vector,
and $\hat{n}$ the direction of the polarization vector $\mathbf{n}$
of the ionizing light. Note that the summation over components $\{l,m,\mu\}$
is coherent, and hence phase sensitive.
\item $Y_{lm}^{*}(\theta_{\hat{k}},\phi_{\hat{k}})$ is a spherical harmonic.
\item $D_{\mu,\mu_{0}}^{1}(R_{\hat{n}})$ is a Wigner rotation matrix element,
with a set of Euler angles $R_{\hat{n}}=(\phi_{\hat{n}},\theta_{\hat{n}},\chi_{\hat{n}})$,
which rotates/projects the polarization into the MF .
\item $I_{\mu_{0}}(\theta_{\hat{k}},\phi_{\hat{k}},\theta_{\hat{n}},\phi_{\hat{n}})$
is the final (observable) MFPAD, for a polarization $\mu_{0}$ and
summed over all symmetry components of the initial and final states,
$\mu_{i}$ and $\mu_{f}$. Note that this sum can be expressed as
an incoherent summation, since these components are (by definition)
orthogonal.
\item $g_{p_{i}}$ is the degeneracy of the state $p_{i}$.
\end{itemize}
The dipole matrix element of eqn. \ref{eq:I} - the radial part of
the dipole matrix element - effectively defines the final state amplitude
and phase. Hence, is equivalent to the general form of eqn. \ref{eq:dipole-int},
but here expanded in terms of symmetries of the light-matter system.

"""

# Imports
import numpy as np
import pandas as pd

# Special functions
from scipy.special import sph_harm
import spherical_functions as sf
import quaternion

def mfpad(dataIn, thres = 1e-2, inds = {'ip':1,'it':1}, res = 50, R = None, p = 0):
    """
    Inputs
    ------
    ind : dictionary, optional.
        Used for sub-selection of matrix elements from Xarrays.
        Default set for length gauage, single it component only, inds = {'ip':1,'it':'1'}.

    R : list of Euler angles or quaternions, optional.
        Define LF > MF polarization geometry/rotations.
        For default case (R = None), 3 geometries are calculated, corresponding to z-pol, x-pol and y-pol cases.
        Defined by Euler angles (p,t,c) = [0 0 0] for z-pol, [0 pi/2 0] for x-pol, [pi/2 pi/2 0] for y-pol.

    p : int, optional.
        Defines LF polarization state, p = -1...1, default p = 0 (linearly pol light along z-axis).
        TODO: add summation over p for multiple pol states in LF.

    """

    # Define reduced data from selection over all data
    daRed = matEleSelector(dataIn, thres = 1e-2, inds = inds)

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

        # Eugler angles for rotation of LF->MF, set as [0 0 0] for z-pol, [0 pi/2 0] for x-pol, [pi/2 pi/2 0] for y-pol
        pRot = [0, 0, np.pi/2]
        tRot = [0, np.pi/2, np.pi/2]
        cRot = [0, 0, 0]

        # Convert to quaternions
        R =  quaternion.from_euler_angles(pRot, tRot, cRot)


    #**************** Calculate MFPADs

    T = []

    # Loop over pol geoms R
    for Rcalc in R:
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
            daTemp = daRed.sel({'mu':mu}) * sf.Wigner_D_element(Rcalc, 1, mu, 0).conj()
            YlmXre = YlmX.reindex_like(daTemp)
            T.append(YlmXre.conj() * daTemp)  # Output full (l,m,mu) expansion



    return T
