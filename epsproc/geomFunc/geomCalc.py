"""
ePSproc geometric terms/functions

Collection of codes for geometric functions and tensors.

26/02/20    v1  Initial implementation.

"""

# Imports
import numpy as np
import pandas as pd
import xarray as xr

# Fast calc functions
from epsproc.geomFunc.w3jVecMethods import w3jguVecCPU, w3jprange
# from epsproc.geomFunc import w3jVecMethods

# Other geom functions
from epsproc.sphCalc import setPolGeoms, wDcalc

# Util funcs.
from epsproc.geomFunc.geomUtils import genllL, selQNsRow, genllpMatE

# Optional imports
try:
    import sparse
except ImportError as e:
    if e.msg != "No module named 'sparse'":
        raise
    print('* sparse not found, sparse matrix forms not available. ')



# NOW SET in ep.util
# def genllL(Lmin = 0, Lmax = 10, mFlag = True):
#     """
#     Generate quantum numbers for angular momentum contractions (l, lp, L)
#
#     Parameters
#     ----------
#     Lmin, Lmax : int, optional, default 0, 10
#         Integer values for Lmin and Lmax respectively.
#
#     mFlag : bool, optional, default = True
#         m, mp take all values -l...+l if mFlag=True, or =0 only if mFlag=False
#
#     Returns
#     -------
#     QNs : 2D np.array
#         Values take all allowed combinations ['l','lp','L','m','mp','M'] up to l=lp=Lmax, one set per row.
#
#     Examples
#     ---------
#     # Calculate up to Lmax = 2
#     >>> QNs = genllL(Lmax=2)
#     # Use with w3jTable function to calculate Wigner 3j terms
#     >>> w3j = w3jTable(QNs = QNs)
#
#
#     """
#
#     # Set QNs for calculation
#     QNs = []
#     for l in np.arange(Lmin, Lmax+1):
#         for lp in np.arange(Lmin, Lmax+1):
#
#             if mFlag:
#                 mMax = l
#                 mpMax = lp
#             else:
#                 mMax = 0
#                 mpMax = 0
#
#             for m in np.arange(-mMax, mMax+1):
#                 for mp in np.arange(-mpMax, mpMax+1):
#                     for L in np.arange(np.abs(l-lp), l+lp+1):
#                         M = -(m+mp)
#                         QNs.append([l, lp, L, m, mp, M])
#
#     return np.array(QNs)
#

def setPhaseConventions(phaseConvention = 'S', typeList = False):
    """
    Set phase convention/choices for geometric functions.

    20/03/20 - first attempt. Aim to centralise all phase choices here to keep things clean and easy to debug/change.

    Set as dictionary for each term, to be appended to results Xarray.


    Parameters
    ----------

    phaseConvention : optional, str, default = 'S'
        Set phase conventions:
        - 'S' : Standard derivation.
        - 'R' : Reduced form geometric tensor derivation.
        - 'E' : ePolyScat, may have additional changes in numerics, e.g. conjugate Wigner D.

    typeList : optional, bool, default = False
        If true, return list of supported options instead of list of phase choices.

    """

    # Return just typeList if option set - this also defines master list.
    if typeList:
        # Supported types
        typeList = ['S', 'R', 'E']
        return typeList


    # Set master dict to hold choices.
    phaseCons = {'phaseConvention':phaseConvention}

    #*** For EPR tensor
    EPRcons = {}
    if phaseConvention == 'S':
        EPRcons['Rphase'] = True        # Apply (-1)^R phase?
        EPRcons['negR'] = True          # Use -R or +R in 3j?
        EPRcons['negRlabel'] = False    # Use -R or +R in QN labels? (Will affect Xarray-based calculations.)
        EPRcons['negRcoordSwap'] = False    # Swap -R and +R in QN coords? (Will affect Xarray-based calculations.)

    elif phaseConvention == 'R':
        EPRcons['Rphase'] = True       # Apply (-1)^R phase?
        EPRcons['negR'] = True          # Use -R or +R in 3j?
        EPRcons['negRlabel'] = False    # Use -R or +R in QN labels? (Will affect Xarray-based calculations.)
        EPRcons['negRcoordSwap'] = False    # Swap -R and +R in QN coords? (Will affect Xarray-based calculations.)

    elif phaseConvention == 'E':
        EPRcons['Rphase'] = True        # Apply (-1)^R phase?
        EPRcons['negR'] = True          # Use -R or +R in 3j?
        EPRcons['negRlabel'] = False    # Use -R or +R in QN labels? (Will affect Xarray-based calculations.)
        EPRcons['negRcoordSwap'] = False    # Swap -R and +R in QN coords? (Will affect Xarray-based calculations.)

    phaseCons['EPR'] = EPRcons

    #*** For Lambda term (as set by MFproj())
    lambdaCons = {}
    if phaseConvention == 'S':
        lambdaCons['negMup'] = True     # Use -mup term in 3j?
        lambdaCons['phaseNegR'] = True  # Set for (-Rp, -R) phase convention (in Wigner D term)?
        lambdaCons['conjFlag'] = False  # Set for conjuate Wigner D?
        lambdaCons['RpPhase'] = True    # Apply (-1)^Rp phase term?

    elif phaseConvention == 'R':
        lambdaCons['negMup'] = True     # Use -mup term in 3j?
        lambdaCons['phaseNegR'] = True  # Set for (-Rp, -R) phase convention (in Wigner D term)?
        lambdaCons['conjFlag'] = False  # Set for conjuate Wigner D?
        lambdaCons['RpPhase'] = True    # Apply (-1)^Rp phase term?

    elif phaseConvention == 'E':
        lambdaCons['negMup'] = True     # Use -mup term in 3j?
        lambdaCons['phaseNegR'] = True  # Set for (-Rp, -R) phase convention (in Wigner D term)?
        lambdaCons['conjFlag'] = True  # Set for conjuate Wigner D?
        lambdaCons['RpPhase'] = True    # Apply (-1)^Rp phase term?

    phaseCons['lambdaCons'] = lambdaCons


    #*** For Beta term (as set by betaTerm())
    betaCons = {}
    if phaseConvention == 'S':
        betaCons['negM'] = False       # Use -M term in 3j?
        betaCons['mPhase'] = True     # Apply (-1)^m phase term?

    elif phaseConvention == 'R':
        betaCons['negM'] = False       # Use -M term in 3j?
        betaCons['mPhase'] = True     # Apply (-1)^m phase term?

    elif phaseConvention == 'E':
        betaCons['negM'] = True       # Use -M term in 3j?
        betaCons['mPhase'] = True     # Apply (-1)^m phase term?

    phaseCons['betaCons'] = betaCons


    #*** For MFPAD product case, as calculated in mfblmXprod()
    mfblmCons = {}
    if phaseConvention == 'S':
        mfblmCons['negRcoordSwap'] = True       # Swap -R and +R in EPR Xarray coords?

        mfblmCons['negMcoordSwap'] = True       # Swap +/-M coords.
        mfblmCons['Mphase'] = True              # Apply (-1)^M phase term.

        mfblmCons['negmCoordSwap'] = True       # Swap +/-m coords.
        mfblmCons['mPhase'] = True              # Apply (-1)^m phase term.

        mfblmCons['mupPhase'] = True            # Apply (-1)^(mup - p) phase term.

        mfblmCons['BLMmPhase'] = False          # TESTING ONLY - switch signs (m, M) terms before 3j calcs.

    if phaseConvention == 'R':
        mfblmCons['negRcoordSwap'] = True       # Swap -R and +R in EPR Xarray coords?

        mfblmCons['negMcoordSwap'] = True       # Swap +/-M coords.
        mfblmCons['Mphase'] = False              # Apply (-1)^M phase term.

        mfblmCons['negmCoordSwap'] = True       # Swap +/-m coords.
        mfblmCons['mPhase'] = True              # Apply (-1)^m phase term.

        mfblmCons['mupPhase'] = True            # Apply (-1)^(mup - p) phase term.

        mfblmCons['BLMmPhase'] = False          # TESTING ONLY - switch signs (m, M) terms before 3j calcs.

    if phaseConvention == 'E':
        mfblmCons['negRcoordSwap'] = True       # Swap -R and +R in EPR Xarray coords?

        mfblmCons['negMcoordSwap'] = False       # Swap +/-M coords.
        mfblmCons['Mphase'] = False              # Apply (-1)^M phase term.

        mfblmCons['negmCoordSwap'] = False       # Swap +/-m coords.
        mfblmCons['mPhase'] = False              # Apply (-1)^m phase term.

        mfblmCons['mupPhase'] = True            # Apply (-1)^(mup - p) phase term.

        mfblmCons['BLMmPhase'] = False          # TESTING ONLY - switch signs (m, M) terms before 3j calcs.

    phaseCons['betaCons'] = betaCons

    return phaseCons


def remapllpL(dataIn, QNs, form = 'dict', method = 'sel', dlist = ['l','lp','L','m','mp','M'], verbose = 0):
    """
    Remap Wigner 3j table, with QNs (l,lp,L,m,mp,M) > tensor forms.

    Tensors are stored by (l,lp,L) triples, with corresponding arrays [m,mp,M], as a dictionary, or Xarray dataarray or dataset.

    Parameters
    ----------

    dataIn : np.array
        Array of data values corresponding to rows (coords) in QNs.

    QNs : np.array
        List of QNs [l, lp, L, m, mp, M] to compute 3j terms for.
        If not supplied all allowed terms for {Lmin, Lmax} will be set.
        (If supplied, values for Lmin, Lmax and mFlag are not used.)

    form : str, optional, default = 'dict'
        Type of output structure.
        - 'dict' : dictionary with keys (l,lp,L), coordinate tables
        - '3d' : dictionary with keys (l,lp,L), 3D arrays indexed by [l+m, lp+mp, L+M]; this case also sets (0,0,0) term as 'w3j0'.
        - 'xdaLM' : Xarray dataarray, with stacked dims ['lSet','mSet']
        - 'xds' : Xarray dataset, with one array per (l,lp,L)
        - 'xdalist' : List of Xarray dataarrays, one per (l,lp,L)

    method : str, optional, default = 'sel'
        Method for selection from input array.
        - 'sel' : use full selection routine, :py:func:`ep.selQNsRow()`, should be most robust.
        - 'n' : sort based on np.unique reindexing. Should be faster, but possibly not robust...?

    dlist : list, optional, default = ['l','lp','L','m','mp','M']
        Full list of dimension names to use/set.

    Returns
    -------
    w3j : Xarray, dictionary
        Wigner3j(l,lp,L,m,mp,M) values corresponding to rows (coords) in QNs, type according to input form.
        Data structure sorted by (l,lp,L) triples.

    """

    # May want to add type checking here...?  QNs should be ints.

    # Find unique (l,lp,L) and corresponding rows
    uniquellpL, setllpL, indllpL = np.unique(QNs[:,:3], axis = 0, return_inverse = True, return_index = True)
    uniquellpL = uniquellpL.astype('int')
    llpLKeys = tuple(map(tuple, uniquellpL))  # Map values to tuple for dict keys

    # Set structures
    w3jDict = {}
    w3jXlist = []
    w3jXDict = {}

    # Loop over (l,lp,L) sets and resturcture data.
    for n, lset in enumerate(uniquellpL):
        if method == 'n':
            w3jDict[llpLKeys[n]] = {'lSet':lset, 'mTable':QNs[indllpL==n, 3:], 'w3jm':dataIn[indllpL==n]}  # Map each unique item back to original list, based on n==set of llpL values.
                                                            # May not be robust?

        # Alternatively, reindex by each unique l,lp,L
        if method == 'sel':
            subset, boolMask = selQNsRow(QNs, lset, fields = [0,1,2], verbose = 0)

            # Set as 3D [m,mp,M] array
            if form == '3d':
                # QNs = subset[:,0:-1].astype('int')  # QNs as ints for indexing

                Lmax = lset.max().astype('int')  # Current Lmax... or set independtly per dim
                # w3jmmpM = np.zeros([2*(Lmax+1), 2*(Lmax+1), 4*(Lmax+1)])
                w3jmmpM = np.zeros([2*lset[0]+1, 2*lset[1]+1, 2*lset[2]+1])  # NOTE - this will fail if w3jTable(nonzeroFlag = False), since illegal terms will be out of index range.

                w3jmmpM[subset[:,0]+subset[:,3], subset[:,1]+subset[:,4], subset[:,2]+subset[:,5]] = dataIn[boolMask]

                # Also set (m,mp,M) = 0 term by index
                indM0 = (np.ceil(np.asarray(w3jmmpM.shape)/2)-1).astype('int')  # Calculate (m,mp,M)=0 coord

                w3jDict[llpLKeys[n]] = {'lSet':lset, 'w3jm':w3jmmpM, 'w3j0':w3jmmpM[indM0[0],indM0[1],indM0[2]]}



            # Test out Xarray datasets, with one entry per (l,lp,L)
            # This allows for sorting by (l,lp,L), then stacking as desired
            elif form.startswith('x'):
                # QNllpL = pd.MultiIndex.from_arrays(lset, names = dlist[0:3])
                QNllpL = pd.MultiIndex.from_tuples([(lset[0],lset[1],lset[2])], names = dlist[0:3])
                QNmmpM = pd.MultiIndex.from_arrays(subset[:,3:].T, names = dlist[3:6])
                # w3jXsub = xr.DataArray(subset[:,-1], coords={'ms':QNmmpM}, dims = ['ms'])
                w3jXsub = xr.DataArray(dataIn[boolMask], coords={'mSet':QNmmpM}, dims = ['mSet'])
                w3jXsub = w3jXsub.expand_dims({'lSet':QNllpL})

                w3jXsub.name = llpLKeys[n]
                w3jXsub.attrs['dataType'] = 'Wigner3j'

                # print(w3jXsub)

                if form.startswith('xda'):
                    w3jXlist.append(w3jXsub)  # List form, can be stacked later
                    # print(w3jXlist)
                else:
                    w3jXDict[llpLKeys[n]] = w3jXsub  # Dict form, can use to construct dataset directly


            # Set as basic coord table
            else:
                w3jDict[llpLKeys[n]] = {'lSet':lset, 'mTable':subset[:,3:], 'w3jm':dataIn[boolMask]}


    # Return values as required
    if (form == 'dict') or (form == '3d'):
        return w3jDict

    elif form == 'xdaLM':
        # Stack list of dataarrays
        daOut = xr.combine_nested(w3jXlist, concat_dim=['lSet'])  # To combine list
        daOut.name = 'w3jStacked'
        return daOut

    elif form == 'xdalist':
        return w3jXlist

    elif form == 'xds':
        # Test dataset - need to pass dataarrays as dict (?) - might be a neater way here...?
        return xr.Dataset(w3jXDict)

    else:
        print(f"Return type {form} not recognised.")
        return None







# Tabulate Wigner 3j terms for a given problem/set of QNs
def w3jTable(Lmin = 0, Lmax = 10, QNs = None, mFlag = True, nonzeroFlag = False, form = '2d', dlist = ['l','lp','L','m','mp','M'], backend = 'par', verbose = 0):
    """
    Calculate/tabulate all wigner 3j terms for a given problem/set of QNs.

.. math::
    \begin{equation}
    \begin{array}{ccc}
    l & l' & L\\
    m & m' & M
    \end{array}
    \end{equation}

    Where l, l' take values Lmin...Lmax (default 0...10).
    :math:`\l-lp\<=L<=l+lp`
    m, mp take values -l...+l if mFlag=True, or =0 only if mFlag=False

    Parameters
    ----------
    Lmin, Lmax : int, optional, default 0, 10
        Integer values for Lmin and Lmax respectively.

    QNs : np.array, optional, default = None
        List of QNs [l, lp, L, m, mp, M] to compute 3j terms for.
        If not supplied all allowed terms for {Lmin, Lmax} will be set.
        (If supplied, values for Lmin, Lmax and mFlag are not used.)

    mFlag : bool, optional, default = True
        m, mp take all values -l...+l if mFlag=True, or =0 only if mFlag=False

    nonzeroFlag : bool, optional, default = False
        Drop null terms before returning values if true.

    form : string, optional, default = '2d'
        Defines return format. Options are:
            - 2d, return 2D np.array, rows [l, lp, L, m, mp, M, 3j]
            - xarray, return xarray
                This is nice for easy selection/indexing, but may be problematic for large Lmax if unstacked (essentailly similar to nd case).
            - nd, return ND np.array, dims indexed as [l, lp, L, l+m, lp+mp, L+M], with values 3j.
                This is suitable for direct indexing, but will have a lot of zero entries and may be large.
            - ndsparse, return ND sparse array, dims indexed as [l, lp, L, l+m, lp+mp, L+M], with values 3j.

        Additional options are set via :py:func:`remapllpL()`. This additionally sorts values by (l,lp,L) triples, which is useful in some cases.
            - 'dict' : dictionary with keys (l,lp,L), coordinate tables
            - '3d' : dictionary with keys (l,lp,L), 3D arrays indexed by [l+m, lp+mp, L+M]; this case also sets (0,0,0) term as 'w3j0'.
            - 'xdaLM' : Xarray dataarray, with stacked dims ['lSet','mSet']
            - 'xds' : Xarray dataset, with one array per (l,lp,L)
            - 'xdalist' : List of Xarray dataarrays, one per (l,lp,L)

    dlist : list of labels, optional, default ['l','lp','L','m','mp','M']
        Used to label array for Xarray output case.

    backend : str, optional, default = 'par'
        See Implementation note below.


    Returns
    -------
    w3j : np.array, Xarray, dictionary
        Wigner3j(l,lp,L,m,mp,M) values corresponding to rows (coords) in QNs, type according to input form.


    Implementation
    --------------
    Currently set to run:
        - 'vec': :py:func:`w3jguVecCPU()`, which uses sf.Wigner3j on the back-end, with additional vectorisation over supplied QNs via Numba's @guvectorize.
        - 'par': :py:func:`w3jprange()`, which uses sf.Wigner3j on the back-end, with parallelization over QNs via Numba's @njit with a prange loop.


    TODO
    -----
    - Move dlist to a utility function.



    """

    # Set QNs, if not supplied.
    if QNs is None:
        QNs = genllL(Lmin = Lmin, Lmax = Lmax, mFlag = mFlag)


    # Calculate 3js using Numba GUvec version
    # w3j_QNs = np.zeros(QNs.shape[0])
    # w3jguVecCPU(QNs, w3j_QNs)

    # Calculate with numba prange version (fastest in testing Feb 2020)
    w3j_QNs = w3jprange(QNs)

    if nonzeroFlag:
        # Drop zero terms
        nzInd = np.nonzero(w3j_QNs)
        w3j_QNs = w3j_QNs[nzInd]
        QNs = QNs[nzInd]

    if form == '2d':

        return np.c_[QNs, w3j_QNs]


    elif form == 'xarray':
        # Multindex case - this retains same size as original, since coords are stacked

        QNs = pd.MultiIndex.from_arrays(QNs.real.T.astype('int8'), names = dlist)
        w3jX = xr.DataArray(w3j_QNs, coords={'QN':QNs}, dims = ['QN'])
        w3jX.attrs['dataType'] = 'Wigner3j'
        return w3jX

    elif form == 'xarray2':
    #     # Multindex case - try statcking by (l,lp,L) and (m,mp,M) sets

        # Currently not working, need to fix dims! .expand_dims might be the way to go?
    #     QNllpL = pd.MultiIndex.from_arrays(QNs[:, 0:3].T.astype('int8'), names = dlist[0:3])
    #     QNmmpM = pd.MultiIndex.from_arrays(QNs[:, 3:6].T.astype('int8'), names = dlist[3:6])
    #     w3jX2 = xr.DataArray(w3j_QNs, coords={'ls':QNllpL, 'ms':QNmmpM}, dims = ['ls','ms'])

        # Array by stack/unstack - works, but will be slow for large arrays however
        QNs = pd.MultiIndex.from_arrays(QNs.real.T.astype('int8'), names = dlist)
        w3jX = xr.DataArray(w3j_QNs, coords={'QN':QNs}, dims = ['QN'])
        w3jX = w3jX.unstack('QN').stack(ls = ['l','lp','L'], ms = ['m','mp','M']).dropna('ms', how = 'all').dropna('ls', how='all')
        w3jX.attrs['dataType'] = 'Wigner3j'
        # Array by resort then stack - avoids large array mem issues, see test dictionary code for method.

        return w3jX


    # Create ND array for direct indexing - this may get large however
    elif form == 'nd':
        w3jND = np.zeros([Lmax+1, Lmax+1, 2*Lmax+1, 2*Lmax+1, 2*Lmax+1, 4*Lmax+1])  # Define full dims

        # Assign values
        # for row in w3j_QNs
        w3jND[QNs[:,0], QNs[:,1], QNs[:,2], QNs[:,0]+QNs[:,3], QNs[:,1]+QNs[:,4], QNs[:,2]+QNs[:,5]] = w3j_QNs
        return w3jND

    # Set to sparse matrix form as defined by Sparse library.
    # https://sparse.pydata.org/
    elif form == 'ndsparse':
        # print('Not implemented')
        w3js = sparse.COO(QNs.T, w3j_QNs)  # Set matrix from ND cood list + data array.
        return w3js

    else:
        # Try additional dict-based outputs in remapllpL
        # Will return none if form is not supported.
        # Note this may involve additional sorting, so be slow for large arrays.
        return remapllpL(w3j_QNs, QNs, form = form, dlist = dlist)





def EPR(QNs = None, p = None, ep = None, nonzeroFlag = True, form = '2d', dlist = ['l', 'lp', 'P', 'p', 'R-p', 'R'], phaseConvention = 'S'):
    """Define polarization tensor (LF) for 1-photon case.

    Define field terms (from QM book, corrected vs. original S\&U version
    - see ``beta-general-forms\_rewrite\_290917.lyx``):

.. math::
    \begin{equation}
    E_{PR}(\hat{e})=[e\otimes e^{*}]_{R}^{P}=[P]^{\frac{1}{2}}\sum_{p}(-1)^{R}\left(\begin{array}{ccc}
    1 & 1 & P\\
    p & R-p & -R
    \end{array}\right)e_{p}e_{R-p}^{*}\label{eq:EPR-defn-1}
    \end{equation}


    Parameters
    ----------
    QNs : np.array, optional, default = None
        List of QNs [l, lp, P, m, mp, R] to compute 3j terms for.
        If not supplied all allowed terms for the one-photon case, l=lp=1, will be set.

    p : list or array, optional, default = None
        Specify polarization terms p.
        If set to None, all allowed values for the one-photon case will be set, p=[-1,0,1]

    ep : list or array, optional, default = None
        Relative strengths for the fields ep.
        If set to None, all terms will be set to unity, ep = 1

    nonzeroFlag : bool, optional, default = True
        Drop null terms before returning values if true.

    form : string, optional, default = '2d'
        For options see :py:func:`ep.w3jTable()`

    phaseConvention : optional, str, default = 'S'
        Set phase conventions:
        - 'S' : Standard derivation.
        - 'R' : Reduced form geometric tensor derivation.
        - 'E' : ePolyScat, may have additional changes in numerics, e.g. conjugate Wigner D.
        See :py:func:`setPhaseConventions` for more details.

    Examples
    ---------
    # Generate full EPR list with defaults
    >>> EPRtable = EPR()

    # Return as Xarray
    >>> EPRtable = EPR(form = 'xarray')

    Note
    ----
    Currently not handling ep correctly!  Should implement as passed Xarray for correct (p,p') assignment.

    """

    # Set phase conventions
    phaseCons = setPhaseConventions(phaseConvention = phaseConvention)

    if phaseCons['EPR']['negRlabel']:
        dlist[-1] = '-R'    # Set -R in QN labels.  Note this assumes dlist[-1] = 'R'

    # Set dim labels for reference/Xarray case - now passed
    # dlist = ['l', 'lp', 'P', 'p', 'R-p', 'R']

    # If p is not passed, set for all allowed values
    if p is None:
        p = [-1,0,1]

    # If no field strength components are defined, set to unity
    if ep is None:
        ep = np.ones(len(p))

    # If no QNs are passed, set for all possible terms
    if QNs is None:
        QNs = []

    # Vanilla version (arb QNs)
    #    for l in np.arange(0, 2):
    #        for lp in np.arange(0, 2):
    #            for m in np.arange(-l, l+1):
    #                for mp in np.arange(-lp, lp+1):
    #                    for L in np.arange(0, l+lp+1):
    #                        M = -(m+mp)
    #                        QNs.append([l, lp, L, m, mp, M])

    # For EPR specific case...
#        for l in np.arange(0, 2):
#            for lp in np.arange(0, 2):
#                for m in np.arange(-l, l+1):
#                    for mp in np.arange(-lp, lp+1):
#                        for R in np.arange(-(l+lp)-1, l+lp+2):
#                            for L in np.arange(0, l+lp+1):
#                                QNs.append([l, lp, L, m, R-mp, -R])
#
    # For case with p preset.
        l = 1
        lp = 1
        for m in p:
            mp = m
            #for R in np.arange(-(l+lp)-1, l+lp+2):
            #    for P in np.arange(0, l+lp+1):
            for P in np.arange(0, l+lp+1):
                for R in np.arange(-P, P+1):

                    # Set phase choice for 3j term
                    if phaseCons['EPR']['negR']:
                        QNs.append([l, lp, P, m, R-mp, -R])
                    else:
                        QNs.append([l, lp, P, m, R-mp, R])

        QNs = np.array(QNs)


    #********** Calculate EPR terms from QN list with Numba guvec function
    # Version with w3jTable()
    EPRtable = w3jTable(QNs = QNs, nonzeroFlag = nonzeroFlag, form = form, dlist = dlist)

    # Phase & degen terms
    if form == '2d':
        Rphase = np.power(-1, np.abs(EPRtable[:,5]))
        Pdegen = np.sqrt(2*EPRtable[:,2] + 1)
        EPRtable[:,-1] = EPRtable[:,-1]*Rphase*Pdegen

    elif form == 'xarray':
        Pdegen = np.sqrt(2*EPRtable.P + 1)

        # Set phase choice for 3j term
        if phaseCons['EPR']['Rphase']:
            Rphase = np.power(-1, np.abs(EPRtable.R))
            EPRtable *= Rphase*Pdegen
        else:
            EPRtable *= Pdegen

        # Switch coord sign?
        if phaseCons['EPR']['negRcoordSwap']:
            temp = EPRtable.unstack()   # Need to unstack to change MultiIndex coords?  Might be a cleaner way?
            temp[dlist[-1]] *= -1       # Label from dlist to allow for +/-R label here.
            EPRtable = temp.stack({'QN':dlist}).dropna(dim = 'QN',how = 'all')
            # NOTE: Without dropna here dims grow! Default settings have 18 elements, but end up with 135 and lots of NaNs.

        EPRtable.attrs['dataType'] = 'EPR'
        EPRtable.attrs['phaseCons'] = phaseCons

    return EPRtable


# BetaTerm for BLM tensor coupling.
def betaTerm(QNs = None, Lmin = 0, Lmax = 10, nonzeroFlag = True, form = '2d', dlist = ['l', 'lp', 'L', 'm', 'mp', 'M'], phaseConvention = 'S'):
    """Define BLM coupling tensor

    Define field terms (from QM book, corrected vs. original S\&U version
    - see ``beta-general-forms\_rewrite\_290917.lyx``):

.. math::
    \begin{equation}
    B_{L,M}=(-1)^{m}\left(\frac{(2l+1)(2l'+1)(2L+1)}{4\pi}\right)^{1/2}\left(\begin{array}{ccc}
    l & l' & L\\
    0 & 0 & 0
    \end{array}\right)\left(\begin{array}{ccc}
    l & l' & L\\
    -m & m' & M
    \end{array}\right)
    \end{equation}

    Parameters
    ----------
    QNs : np.array, optional, default = None
        List of QNs [l, lp, L, m, mp, M] to compute 3j terms for.
        If not supplied all allowed terms for {Lmin, Lmax} will be set.
        (If supplied, values for Lmin, Lmax and mFlag are not used.)

    Lmin, Lmax : int, optional, default 0, 10
        Integer values for Lmin and Lmax respectively.

    mFlag : bool, optional, default = True
        m, mp take all values -l...+l if mFlag=True, or =0 only if mFlag=False

    nonzeroFlag : bool, optional, default = True
        Drop null terms before returning values if true.

    form : string, optional, default = '2d'
        For options see :py:func:`ep.w3jTable()`

    dlist : list of labels, optional, default ['l','lp','L','m','mp','M']
        Used to label array for Xarray output case.

    phaseConvention : optional, str, default = 'S'
        Set phase conventions:
        - 'S' : Standard derivation.
        - 'R' : Reduced form geometric tensor derivation.
        - 'E' : ePolyScat, may have additional changes in numerics, e.g. conjugate Wigner D.
        See :py:func:`setPhaseConventions` for more details.

    Examples
    ---------

    >>> Lmax = 2
    >>> BLMtable = betaTerm(Lmax = Lmax, form = 'xds')
    >>> BLMtable = betaTerm(Lmax = Lmax, form = 'xdaLM')

    """

    # Define phase conventions for different forms of the term
    phaseCons = setPhaseConventions(phaseConvention = phaseConvention)

    # Set -M for 3j term if required - only if QNs passed however.
    if phaseCons['betaCons']['negM'] and (QNs is not None):
        QNs[:,-1] *= -1

    # Tabulate 3j terms
    BLMtable = w3jTable(QNs = QNs, Lmin = Lmin, Lmax = Lmax, nonzeroFlag = nonzeroFlag, form = form, dlist = dlist)

    # Phase & degen terms
    # Currently only implemented for some cases.

#     if form == '2d':
#         Rphase = np.power(-1, np.abs(EPRtable[:,5]))
#         Pdegen = np.sqrt(2*EPRtable[:,2] + 1)
#         EPRtable[:,-1] = EPRtable[:,-1]*Rphase*Pdegen

#     elif form == 'xarray':
#         Rphase = np.power(-1, np.abs(EPRtable.R))
#         Pdegen = np.sqrt(2*EPRtable.P + 1)
#         EPRtable *= Rphase*Pdegen

    if (form == 'xdaLM') or (form == 'xds'):

        # 3j product term
        try:
            # 3j product term
            # BLMtable *= mPhase*np.sqrt(degen)*BLMtable.sel(m=0,mp=0,M=0)  # Accidental correlations in results...?
            BLMtable = BLMtable*BLMtable.sel(m=0,mp=0,M=0).drop('mSet').squeeze()  # NOTE - drop dims here to prevent correlated product.

        # If (0,0,0) terms are not already calculated, do so.
        # Might be cleaner just to use this in all cases?
        except KeyError:
            thrj0 = ep.geomFunc.w3jTable(QNs = genllpMatE(matE, mFlag = False), nonzeroFlag = True, form = form, dlist = ['l', 'lp', 'L', 'm', 'mp', 'M'])
            BLMtable = BLMtable*thrj0.drop('mSet').squeeze()

        mPhase = np.power(-1, np.abs(BLMtable.m))
        degen = (2*BLMtable.l+1)*(2*BLMtable.lp+1)*((2*BLMtable.L+1))/(4*np.pi)

        if phaseCons['betaCons']['mPhase']:
            BLMtable *= mPhase*np.sqrt(degen)
        else:
            BLMtable *= np.sqrt(degen)

        BLMtable.attrs['dataType'] = 'betaTerm'
        BLMtable.attrs['phaseCons'] = phaseCons

    else:
        print(f"Form {form} not implemented.")
        BLMtable = None

    return BLMtable


# Define lambdaTerm, MF projection.
def MFproj(QNs = None, RX = None, nonzeroFlag = True, form = '2d', dlist = ['l', 'lp', 'P', 'mu', 'mup', 'Rp', 'R'], phaseConvention = 'S'):
    """
    Define MF projection term, :math:`\Lambda_{R',R}(R_{\hat{n}})`:

.. math::
    \begin{equation}
    \Lambda_{R',R}(R_{\hat{n}})=(-1)^{(R')}\left(\begin{array}{ccc}
    1 & 1 & P\\
    \mu & -\mu' & R'
    \end{array}\right)D_{-R',-R}^{P}(R_{\hat{n}})
    \end{equation}


    Then...

.. math::
    \begin{eqnarray}
    \beta_{L,-M}^{\mu_{i},\mu_{f}} & = & \sum_{P,R',R}{\color{red}E_{P-R}(\hat{e};\mu_{0})}\sum_{l,m,\mu}\sum_{l',m',\mu'}(-1)^{(\mu'-\mu_{0})}{\color{red}\Lambda_{R',R}(R_{\hat{n}};\mu,P,R,R')B_{L,-M}(l,l',m,m')}I_{l,m,\mu}^{p_{i}\mu_{i},p_{f}\mu_{f}}(E)I_{l',m',\mu'}^{p_{i}\mu_{i},p_{f}\mu_{f}*}(E)
    \end{eqnarray}

    Parameters
    ----------

    phaseConvention : optional, str, default = 'S'
        Set phase conventions:
        - 'S' : Standard derivation.
        - 'R' : Reduced form geometric tensor derivation.
        - 'E' : ePolyScat, conjugate Wigner D.
        See :py:func:`setPhaseConventions` for more details.

    Notes
    -----
    This is very similar to :math:`E_{PR}` term.

    Examples
    --------
    >>> lTerm, lambdaTable, lambdaD, QNs = MFproj(form = 'xarray')

    """

    # Define phase conventions for different forms of the term
    # phaseNegR: Set for (-Rp, -R) phase convention (otherwise use (+Rp,+R))
    # conjFlag: Set for conjugate Wigner D terms
    # if phaseConvention == 'E':
    #     phaseNegR = True
    #     conjFlag = True
    #     # QNphase = True
    # else:
    #     phaseNegR = True
    #     conjFlag = False
    #     # QNphase = False
    # Set phase conventions
    phaseCons = setPhaseConventions(phaseConvention = phaseConvention)

    # If no QNs are passed, set for all possible terms
    if QNs is None:
        QNs = []

        # Set photon terms
        l = 1
        lp = 1

        # Loop to set all other QNs
        for mu in np.arange(-l, l+1):
            for mup in np.arange(-lp, lp+1):
                #for R in np.arange(-(l+lp)-1, l+lp+2):
                #    for P in np.arange(0, l+lp+1):
                for P in np.arange(0, l+lp+1):
                    for Rp in np.arange(-P, P+1):
                        for R in np.arange(-P, P+1):
                            # QNs.append([l, lp, P, mu, -mup, R, Rp])
                            if phaseCons['lambdaCons']['negMup']:
                                QNs.append([l, lp, P, mu, -mup, R, Rp])
                            else:
                                QNs.append([l, lp, P, mu, mup, R, Rp])

        QNs = np.array(QNs)

    # Set to calculate for default (x,y,z) pol geoms.
    if RX is None:
        RX = setPolGeoms()

    #*********************** Calculate terms
    # 3j term with w3jTable()
    # lambdaTable = w3jTable(QNs = QNs[:,:-1], nonzeroFlag = nonzeroFlag, form = form, dlist = dlist[:-1])  # Pass only dims for 3j term
    if form == '2d':
        nonzeroFlag = False   # For 2d case override to ensure consistent size.

    lambdaTable = w3jTable(QNs = QNs, nonzeroFlag = nonzeroFlag, form = form, dlist = dlist)  # Pass all QNs to keep R label. w3jpRange will select/use QN[:,0:5] only

    # D term with wDcalc(), includes order switch for (R,Rp)
    # Pass RX.data since wDcalc is curently only set for np.array inputs.
    # Set XFlag by form.
    if form.startswith('x'):
        XFlag = True
    else:
        XFlag = False

    # Subselect on QNs and calculate
    QNind = np.array([2, 6, 5])   # Set (P,Rp,R)
    dRed = [dlist[n] for n in QNind]

    QNwD = QNs[:,QNind]
    QNun = np.unique(QNs[:,QNind], axis=0)

    # Set for (-Rp, -R) phase convention
    if phaseCons['lambdaCons']['phaseNegR']:
        QNwD[:,1:] *= -1
        QNun[:,1:] *= -1

    if form == '2d':
        # Cal for all values - in cases with duplicate QNs this may lead to indexing issues later in Xarray case.
        # Should be OK for 2d case however, will provide table matching full QN list.
        # NOTE in 2d case, wDcalc() also outputs R, QNs - skip these since they're already set in this case
        lambdaD, *RQN = wDcalc(QNs = QNwD, R = RX.data, XFlag = XFlag, dlist = dRed, eNames = ['Phi','Theta','Chi'],
                                conjFlag = phaseCons['lambdaCons']['conjFlag'])
        lambdaD = np.asarray(lambdaD)

    elif form.startswith('x'):
        # Calc for unique values only to avoid duplicate coords
        lambdaD = wDcalc(QNs = QNun, R = RX.data, XFlag = XFlag, dlist = dRed, eNames = ['Phi','Theta','Chi'],
                         conjFlag = phaseCons['lambdaCons']['conjFlag'])

        lambdaD['Labels']=('Euler',RX.Labels.values)  # Propagate labels, currently wDcalc only takes RX.data
        lambdaD = lambdaD.swap_dims({'Euler':'Labels'})  # Swap dims to labels.
        lambdaD.attrs['dataType'] = 'WignerD'



    #***************** Multiplications & phase
    if form.startswith('x'):
        lun = lambdaTable.unstack('QN')
        lDun = lambdaD.unstack('QN')

        # Reset phase choices here to allow for correct multiplication of +Rp and D(P,-Rp,-R) terms in Xarray
        if phaseCons['lambdaCons']['phaseNegR']:
            lDun['R'] *= -1
            lDun['Rp'] *= -1

        # Additional phase term (-1)^(-Rp)
        if phaseCons['lambdaCons']['RpPhase']:
            Rpphase = np.power(-1, np.abs(lun.Rp))
            lTerm = Rpphase * lun * lDun
        else:
            lTerm = lun * lDun

        lTerm.attrs['dataType'] = 'LambdaTerm'
        lTerm.attrs['phaseCons'] = phaseCons

    # Array multiplication case
    elif form == '2d':

        # Additional phase term (-1)^(-Rp)
        if phaseCons['lambdaCons']['RpPhase']:
            Rpphase = np.power(-1, np.abs(QNs[:,-1]))
        else:
            Rpphase = np.ones(QNs[:,-1].size)

        # Loop over sets of Euler angles
        lTerm = []
        for eInd in range(lambdaD.shape[1]):
            lTerm.append(Rpphase * lambdaD[:,eInd] * lambdaTable[:,-1])

        lTerm = np.array(lTerm).T

    else:
        print(f'Form {form} not supported.')


    return lTerm, lambdaTable, lambdaD, QNs
