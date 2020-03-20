# -*- coding: utf-8 -*-
"""
ePSproc spherical function calculations.

Collection of functions for calculating Spherical Tensors: Ylm, wignerD etc.

For spherical harmonics, currently using scipy.special.sph_harm

For other functions, using Moble's spherical_functions package
https://github.com/moble/spherical_functions

See tests/Spherical function testing Aug 2019.ipynb

04/12/19        Added `setPolGeoms()` to define frames as Xarray.
                Added `setADMs()` to define ADMs as Xarray
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
import string

try:
    from sympy.physics.quantum.spin import Rotation  # For basic frame rotation code, should update to use sf
except ImportError as e:
    if e.msg != "No module named 'sympy'":
        raise
    print('* Sympy not found, some (legacy) sph functions may not be available. ')



# Master function for setting geometries/frame rotations
def setPolGeoms(eulerAngs = None, quat = None, labels = None, vFlag = 2):
    """
    Generate Xarray containing polarization geometries as Euler angles and corresponding quaternions.

    Define LF > MF polarization geometry/rotations.
    Provide either eulerAngs or quaternions, but not both (supplied quaternions only will be used in this case).

    For default case (eulerAngs = None, quat = None), 3 geometries are calculated,
    corresponding to z-pol, x-pol and y-pol cases.
    Defined by Euler angles:
    (p,t,c) = [0 0 0] for z-pol,
    (p,t,c) = [0 pi/2 0] for x-pol,
    (p,t,c) = [pi/2 pi/2 0] for y-pol.


    Parameters
    ----------
    eulerAngs : list or np.array of Euler angles (p(hi), t(heta), c(hi)), optional.
        List or array [p,t,c...], shape (Nx3).
        List or array including set labels, [label,p,t,c...], shape (Nx4)

    quat : list or np.array of quaternions, optional.

    labels : list of labels, one per set of angles. Optional.
        If not set, states will be labelled numerically.

    vFlag : version of routine to use, optional, default = 2
        Options:
        - 1, use labels as sub-dimensional coord.
        - 2, set labels as non-dimensional coord.

    Returns
    -------
    RX : Xarray of quaternions, with Euler angles as dimensional params.

    To do
    -----
    - Better label handling, as dictionary? With mixed-type array may get issues later.
        (sf.quaternion doesn't seem to have an issue however.)
    - Xarray MultiIndex with mixed types?
        Tested with pd - not supported:
        >>> eulerInd = pd.MultiIndex.from_arrays([eulerAngs[:,0].T, eulerAngs[:,1:].T.astype('float')], names = ['Label','P','T','C'])
        # Gives error:
        #   NotImplementedError: > 1 ndim Categorical are not supported at this time


    Examples
    --------

    >>> # Defaults
    >>> RXdefault = setPolGeoms()
    >>> print(RXdefault)

    >>> # Pass Eulers, no labels
    >>> pRot = [1.666, 0, np.pi/2]
    >>> tRot = [0, np.pi/2, np.pi/2]
    >>> cRot = [-1.5, 0, 0]

    >>> eulerAngs = np.array([pRot, tRot, cRot]).T

    >>> RXePass = setPolGeoms(eulerAngs = eulerAngs)
    >>> print(RXePass)

    >>> # Pass labels separately
    >>> RXePass = setPolGeoms(eulerAngs = eulerAngs, labels = ['1','23','ff'])
    >>> print(RXePass)

    >>> # Pass Eulers with existing labels
    >>> labels = ['A','B','C']
    >>> eulerAngs = np.array([labels, pRot, tRot, cRot]).T
    >>> RXePass = setPolGeoms(eulerAngs = eulerAngs)
    >>> print(RXePass)

    >>> # Pass Quaternions and labels
    >>> RXqPass = setPolGeoms(quat = RXePass, labels = labels)
    >>> print(RXqPass)

    >>> # Pass both - only quaternions will be used in this case, and warning displayed.
    >>> RXqeTest = setPolGeoms(eulerAngs = eulerAngs, quat = RXePass, labels = labels)
    >>> print(RXqeTest)

    """

    # Default case, set (x,y,z) geometries
    if (eulerAngs is None) and (quat is None):
        # As arrays, with labels
        pRot = [0, 0, np.pi/2]
        tRot = [0, np.pi/2, np.pi/2]
        cRot = [0, 0, 0]
        labels = ['z','x','y']
        eulerAngs = np.array([labels, pRot, tRot, cRot]).T   # List form to use later, rows per set of angles

    # Get quaternions from Eulers, if provided or as set above for default case.
    if eulerAngs is not None:
        if type(eulerAngs) is not np.ndarray:
            eulerAngs = np.asarray(eulerAngs)

        if eulerAngs.shape[1] is 3:
            if labels is None:
                # Set labels if missing, alphabetic or numeric
                if eulerAngs.shape[0] < 27:
                    labels = list(string.ascii_uppercase[0:eulerAngs.shape[0]])
                else:
                    labels = np.arange(1,eulerAngs.shape[0]+1)

            eulerAngs = np.c_[labels, eulerAngs]

    # If quaternions are passed, set corresponding Eulers
    if quat is not None:
        eulerFromQuat = quaternion.as_euler_angles(quat) # Set Eulers from quaternions

        if labels is None:
            # Set labels if missing
            labels = np.arange(1,eulerFromQuat.shape[0]+1)

        if eulerAngs is not None:
            print('***Warning: Euler angles and Quaternions passed, using Quaternions only.')

        eulerAngs = np.c_[labels, eulerFromQuat]

    # Otherwise set from Eulers
    else:
        quat = quaternion.from_euler_angles(eulerAngs[:,1:]) # Convert Eulers to quaternions


    #*** Set up Xarray

    if vFlag == 1:
        # v1    keep Labels as subdim.
        #       This works, and allows selection by label, but Euler coords may be string type
        # Set Pandas MultiIndex - note transpose for eulerAngs to (angs,set) order
        eulerInd = pd.MultiIndex.from_arrays(eulerAngs.T, names = ['Label','P','T','C'])
        # Create Xarray
        RX = xr.DataArray(quat, coords={'Euler':eulerInd}, dims='Euler')
        RX.attrs['dataType'] = 'Euler'

    elif vFlag == 2:
        # v2    Labels as non-dim coords.
        #       Doesn't allow selection, but keeps Euler coords as floats in all cases.
        Euler = pd.MultiIndex.from_arrays(eulerAngs[:,1:].T.astype('float'), names = ['P','T','C'])
        RX = xr.DataArray(quat, coords={'Euler':Euler,'Labels':('Euler',eulerAngs[:,0].T)}, dims='Euler')
        RX.attrs['dataType'] = 'Euler'

    else:
        print('***Version not recognized')

    return RX


# Create Xarray from set of ADMs - adapted from existing blmXarray()
def setADMs(ADMs = [0,0,0,1], KQSLabels = None, t = None, addS = False):
    """
    Create Xarray from ADMs, or create default case ADM(K,Q,S) = [0,0,0,1].

    Parameters
    ----------
    ADMs : list or np.array, default = [0,0,0,1]
        Set of ADMs = [K, Q, S, ADM].
        If multiple ADMs are provided per (K,Q,S) index, they are set to the t axis (if provided), or indexed numerically.

    KQSLabels : list or np.array, optional, default = None
        If passed, assume ADMs are unabelled, and use (K,Q,S) indicies provided here.

    t : list or np.array, optional, default = None
        If passed, use for dimension defining ADM sets (usually time).
        Defaults to numerical label if not passed, t = np.arange(0,ADMs.shape[1])

    addS : bool, default = False
        If set, append S = 0 to ADMs.
        This allows for passing of [K,Q,ADM] type values (e.g. for symmetric top case)

    Returns
    -------
    ADMX : Xarray
        ADMs in Xarray format, dims as per :py:func:`epsproc.utils.ADMdimList()`

    Examples
    ---------
    >>> # Default case
    >>> ADMX = setADMs()
    >>> ADMX

    >>> # With full N2 rotational wavepacket ADM set from demo data (ePSproc\data\alignment), where modPath defines root...
    >>> # Load ADMs for N2
    >>> from scipy.io import loadmat
    >>> ADMdataFile = os.path.join(modPath, 'data', 'alignment', 'N2_ADM_VM_290816.mat')
    >>> ADMs = loadmat(ADMdataFile)
    >>> ADMX = setADMs(ADMs = ADMs['ADM'], KQSLabels = ADMs['ADMlist'], addS = True)
    >>> ADMX

    """

    # Check size of passed set of ADMs
    # For ease of manipulation, just change to np.array if necessary!
    if isinstance(ADMs, list):
        ADMs = np.array(ADMs, ndmin = 2)

    # Set lables explicitly if not passed, and resize ADMs
    if KQSLabels is None:
        if addS:
            KQSLabels = ADMs[:,0:2]
            KQSLabels = np.c_[KQSLabels, np.zeros(KQSLabels.shape[0])]
            ADMs = ADMs[:,2:]
        else:
            KQSLabels = ADMs[:,0:3]
            ADMs = ADMs[:,3:]
    else:
        if addS:
            KQSLabels = np.c_[KQSLabels, np.zeros(KQSLabels.shape[0])]  # Add S for labels passed case


    # Set indexing, default to numerical
    if t is None:
        t = np.arange(0,ADMs.shape[1])


    # Set up Xarray
    QNs = pd.MultiIndex.from_arrays(KQSLabels.real.T.astype('int8'), names = ['K','Q','S'])  # Set lables, enforce type
    ADMX = xr.DataArray(ADMs, coords={'ADM':QNs,'t':t}, dims = ['ADM','t'])

    ADMX.attrs['dataType'] = 'ADM'

    return ADMX


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
def wDcalc(Lrange = [0, 1], Nangs = None, eAngs = None, R = None, XFlag = True, QNs = None, dlist = ['lp','mu','mu0'], eNames = ['P','T','C'], conjFlag = False):
    '''
    Calculate set of Wigner D functions D(l,m,mp; R) on a grid.

    Parameters
    ----------
    Lrange : list, optional, default [0, 1]
        Range of L to calculate parameters for.
        If len(Lrange) == 2 assumed to be of form [Lmin, Lmax], otherwise list is used directly.
        For a given l, all (m, mp) combinations are calculated.

    QNs : np.array, optional, default = None
        List of QNs [l,m,mp] to compute Wigner D terms for.
        If supplied, use this instead of Lrange setting.

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

    dlist : list, optional, default ['lp','mu','mu0']
        Labels for Xarray QN dims.
    eNames : list, optional, default ['P','T','C']
        Labels for Xarray Euler dims.

    conjFlag : bool, optional, default = False
        If true, return complex conjuage values.

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

    For testing, see https://epsproc.readthedocs.io/en/latest/tests/Spherical_function_testing_Aug_2019.html

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

    # Set QNs based on Lrange if not passed to function.
    if QNs is None:
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
    else:
        if R is not None:
            eAngs = quaternion.as_euler_angles(R) # Set Eulers from quaternions

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

        if conjFlag:
            wD.append(sf.Wigner_D_element(R, QNs[n,0], QNs[n,1], QNs[n,2]).conj())
        else:
            wD.append(sf.Wigner_D_element(R, QNs[n,0], QNs[n,1], QNs[n,2]))

    # Return values as Xarray or np.arrays
    if XFlag:
        # Put into Xarray
        #TODO: this will currently fail for a single set of QNs.
        QNs = pd.MultiIndex.from_arrays(np.asarray(lmmp).T, names = dlist)
        if (eAngs is not None) and (eAngs.size == 3):  # Ugh, special case for only one set of angles.
            Euler = pd.MultiIndex.from_arrays([[eAngs[0]],[eAngs[1]],[eAngs[2]]], names = eNames)
            wDX = xr.DataArray(np.asarray(wD), coords=[('QN',QNs)])
            wDX = wDX.expand_dims({'Euler':Euler})
        else:
            Euler = pd.MultiIndex.from_arrays(eAngs.T, names = eNames)
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

# 05/12/19 Rewriting with new eAngs and ADM defns... (Xarrays)
def TKQarrayRotX(TKQin, RX, form = 2):
    r"""
    Frame rotation for multipoles $T_{K,Q}$.

    Basic frame rotation code, see https://github.com/phockett/Quantum-Metrology-with-Photoelectrons/blob/master/Alignment/Alignment-1.ipynb for examples.

    Parameters
    ----------
    TKQin : Xarray
        Values defining the initial distribution, [K,Q,TKQ]. Other dimensions will be propagated.

    RX : Xarray defining frame rotations, from :py:func:`epsproc.setPolGeoms()`
        List of Euler angles (theta,phi,chi) and corresponding quaternions defining rotated frame.

    Returns
    -------

    TKQRot : Xarray
        Multipoles $T'_{K,Q}$ in rotated frame, as an np.array [K,Q,TKQ].


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

    Examples
    --------

    >>> vFlag = 2
    >>> RX = ep.setPolGeoms(vFlag = vFlag)  # Package version
    >>> RX
    >>> testADMX = ep.setADMs(ADMs=[[0,0,0,1],[2,0,0,0.5]])
    >>> testADMX
    >>> testADMrot, wDX, wDXre = TKQarrayRotX(testADMX, RX)
    >>> testADMrot
    >>> testADMrot.attrs['dataType'] = 'ADM'
    >>> sph, _ = sphFromBLMPlot(testADMrot, facetDim = 'Euler', plotFlag = True)

    """

    # Check dataType and rename if required
    if TKQin.dataType is 'ADM':
        TKQ = TKQin.copy()
    elif TKQin.dataType is 'BLM':
        TKQ = TKQin.copy().unstack('BLM').rename({'l':'K','m':'Q'}).stack({'ADM':('K','Q')})
    else:
        print('***TKQ dataType not recognized, skipping frame rotation.')
        return None, None, None

    # Test if S is set, and flag for later
    # Better way to get MultiIndexes here?
    if 'S' in TKQ.unstack().dims:
#         incS = TKQ.S.pipe(np.abs).values.max() > 0
        incS = True
    else:
        incS = False

    # If S = 0, apply basic TKQ transformation
#     if not incS:

        #*** Formulate using existing style (looped)
#         # Loop over rotations, extract quaternion value from Xarray (better way to do this...?)
#         # Note this will fail for looping over RX, then taking values - seems to give size=1 array which throws errors... weird...
#         for R in RX.values:
#             # Loop over input K values, for all Q
#             for Kin in ADMX.K:

#                 # Set QNs

#                 # Rotation matrix elements
#                 sf.Wigner_D_element(R, QNs)

        #*** Formulate using existing wDX code, then multiply - should be faster and transparent (?), and allow multiple dims

    # Calculate Wigner Ds
    wDX = wDcalc(Lrange = np.unique(TKQ.K.values), R = RX.values)  # NOTE - alternatively can pass angles as Eulers, but may need type conversion for RX.Euler depending on format, and/or loop over angle sets.

        # Rename coords, use dataType for this
#         dtList = ep.dataTypesList()
#         dtDims = dtList[TKQ.dataType]['dims']  # Best way to get relevant QNs here...?  Maybe need to start labelling these in dataTypesList?

    # Rename for ADMs for testing...
    # ... then mutliply (with existing dims), resort & sum over Q
    if incS:

        # Test cases
        if form == 1:
            wDXre = wDX.unstack('QN').rename({'lp':'K','mu':'Q','mu0':'Qp'}).expand_dims({'S':[0]}).stack({'ADM':('K','Q','S')})  # Restack according to D^l_{mu,mu0} > D^K_{Q,Qp}
                                                                                                                                # This matches Blum 4.41 for CONJ(wD) and conj(TKQ)
             # Here formulated as TKQp* = sum_Q(TKQ* x D^K_{Q,Qp}*)
            TKQrot = (TKQ.conj() * wDXre.conj()).unstack('ADM').sum('Q').rename({'Qp':'Q'}).stack({'ADM':('K','Q','S')}).conj()
            # Gives *no difference* between (x,y) cases? Should be phase rotation?

        if form == 2:  # ******************* THINK THIS is the correct case.

            wDXre = wDX.unstack('QN').rename({'lp':'K','mu':'Qp','mu0':'Q'}).expand_dims({'S':[0]}).stack({'ADM':('K','Q','S')})  # Restack according to D^l_{mu,mu0} > D^K_{Qp,Q}
                                                                                                                             # This matches Zare, eqn. 3.83, for D*xTKQ
            # Here formulated as TKQrot = sum_q(TKq x D^K_{q,Q}*)
            TKQrot = (TKQ * wDXre.conj()).unstack('ADM').sum('Q').rename({'Qp':'Q'}).stack({'ADM':('K','Q','S')})
            # Gives re/im difference between (x,y) cases? Should be phase rotation?

        if form == 3:
            wDXre = wDX.unstack('QN').rename({'lp':'K','mu':'Q','mu0':'Qp'}).expand_dims({'S':[0]}).stack({'ADM':('K','Q','S')})  # Restack according to D^l_{mu,mu0} > D^K_{Q,Qp}
                                                                                                                                # This matches Zare, eqn. 5.8, for sum over Q and REAL wD
            # Here formulated as TKQp = sum_Q(TKQ x D^K_{Q,Qp})
            TKQrot = (TKQ * wDXre).unstack('ADM').sum('Q').rename({'Qp':'Q'}).stack({'ADM':('K','Q','S')})
            # TKQrot = (wDXre * TKQ).unstack('ADM').sum('Q').rename({'Qp':'Q'}).stack({'ADM':('K','Q','S')})
            # Gives *no difference* between (x,y) cases? Should be phase rotation?

    else:
        # wDXre = wDX.unstack('QN').rename({'lp':'K','mu':'Q','mu0':'Qp'}).stack({'ADM':('K','Q')})  # Restack according to D^l_{mu,mu0} > D^K_{Q,Qp}
        # TKQrot = (TKQ * wDXre).unstack('ADM').sum('Q').rename({'Qp':'Q'}).stack({'ADM':('K','Q')})

        # form = 2 case only.
        # wDXre = wDX.unstack('QN').rename({'lp':'K','mu':'Qp','mu0':'Q'}).expand_dims({'S':[0]}).stack({'ADM':('K','Q','S')})  # Restack according to D^l_{mu,mu0} > D^K_{Qp,Q}
        wDXre = wDX.unstack('QN').rename({'lp':'K','mu':'Qp','mu0':'Q'}).stack({'ADM':('K','Q')})
                                                                                                 # This matches Zare, eqn. 3.83, for D*xTKQ
        # Here formulated as TKQrot = sum_q(TKq x D^K_{q,Q}*)
        TKQrot = (TKQ * wDXre.conj()).unstack('ADM').sum('Q').rename({'Qp':'Q'}).stack({'ADM':('K','Q')})

    #***  Mutliply (with existing dims), then resort & sum over Q
    # NOW INCLUDED ABOVE for different test cases


    # Propagate frame labels & attribs
    # TODO: fix Labels propagation - this seems to drop sometimes, dim issue?
    # TKQrot['Labels'] = RX.Labels
    TKQrot['Labels']=('Euler',RX.Labels.values)  # This seems to work...
    TKQrot.attrs = TKQ.attrs

    # For BLM data, rename vars.
    if TKQin.dataType is 'BLM':
        TKQrot = TKQrot.unstack('ADM').rename({'K':'l','Q':'m'}).stack({'BLM':('l','m')})

    return TKQrot, wDX, wDXre
