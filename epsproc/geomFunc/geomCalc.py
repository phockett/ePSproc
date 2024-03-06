"""
ePSproc geometric terms/functions

Collection of codes for geometric functions and tensors.

26/02/20    v1  Initial implementation.
18/07/23    v2  Updated phase conventions following 3D AFPAD testing.
                Fixed use of '*=' in-place for Xarray, no longer robust.
                NOTE: changed for index cases only.

"""

# Imports
import numpy as np
import pandas as pd
import xarray as xr

# Fast calc functions
# from epsproc.geomFunc.w3jVecMethods import Wigner3jQNs, w3jguVecCPU, w3jprange
# 11/10/22 - removed w3jguVecCPU, issues with some Numba versions.
from epsproc.geomFunc.w3jVecMethods import Wigner3jQNs, w3jprange
# from epsproc.geomFunc import w3jVecMethods

# Sympy wrappers
from epsproc.geomFunc.wNjSympyWrapper import w3jSympy

# Other geom functions
from epsproc.sphCalc import setPolGeoms, wDcalc

# Util funcs.
from epsproc.geomFunc.geomUtils import genllL, selQNsRow, genllpMatE, genllLList, genKQStermsFromTensors, genLmLamKQStermsFromTensors

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

#*************************************************************
#************** Set and sort functions (see also geomUtils.py)
#*************************************************************

def setPhaseConventions(phaseConvention = 'S', typeList = False):
    """
    Set phase convention/choices for geometric functions.

    20/03/20 - first attempt. Aim to centralise all phase choices here to keep things clean and easy to debug/change.

    July 2023 updates:
        - Set EPRcons['negRcoordSwap'] = True in all cases. This swaps -R and +R in QN coords (Xarray only), and fixes issues with phase mismatches later.
        - Set phaseCons['afblmCons'] per type. This decouples testing cases.
        -

    Set as dictionary for each term, to be appended to results Xarray.


    Parameters
    ----------

    phaseConvention : optional, str, default = 'S'
        Set phase conventions:
        - 'S' : Standard derivation.
        - 'R' : Reduced form geometric tensor derivation.
        - 'E' : ePolyScat, may have additional changes in numerics, e.g. conjugate Wigner D.
        If a dict of phaseConventions is passed they will simply be returned - this is for transparency/consistency over multiple fns which call setPhaseConventions()... although may be an issue in some cases.

    typeList : optional, bool, default = False
        If true, return list of supported options instead of list of phase choices.


    Note
    -----
    If a dict of phaseConventions is passed they will simply be returned - this is for transparency/consistency over multiple fns which call setPhaseConventions()... although may be an issue in some cases.

    """

    # Return just typeList if option set - this also defines master list.
    if typeList:
        # Supported types
        typeList = ['S', 'R', 'E']
        return typeList

    # If phaseConventions are preset, just return them.
    if type(phaseConvention) is dict:
        return phaseConvention


    # Set master dict to hold choices.
    phaseCons = {'phaseConvention':phaseConvention}

    #**** For generating QNs with genllpMatE()
    # Set here to avoid issues with dropped/missing terms later!
    # Some conventions will be tied to other choices below.
    genMatEcons = {}
    if phaseConvention == 'S':
        genMatEcons['negm'] = False     # Set -m, corresponding to M = -m + mp, otherwise M = -(m+mp)
        # 18/07/23 - testing to match phaseCons['afblmCons']['negM'] = not(phaseCons['genMatEcons']['negm'])  # 01/07/23 - set to not() here, should be correct and tested for 3D alignment case (doesn't affect 1D results)
        # If False, get inconsistencies in Delta term output (no M!=0) for test cases.
        # genMatEcons['negm'] = True     # Set -m, corresponding to M = -m + mp, otherwise M = -(m+mp)
        # genMatEcons['negM'] = False     # Set -M

    elif phaseConvention == 'R':
        genMatEcons['negm'] = False     # Set -m, corresponding to M = -m + mp, otherwise M = -(m+mp)
        # genMatEcons['negM'] = False     # Set -M

    elif phaseConvention == 'E':
        genMatEcons['negm'] = False     # Set -m, corresponding to M = -m + mp (normal convention), otherwise M = -(m+mp)
                                        # Note this is correlated with M switch terms later, incorrect settings may remove m or m' non-zero terms!

        # genMatEcons['negM'] = False     # Set -M

    phaseCons['genMatEcons'] = genMatEcons

    #*** For EPR tensor
    EPRcons = {}
    if phaseConvention == 'S':
        EPRcons['Rphase'] = True        # Apply (-1)^R phase?
        EPRcons['negR'] = True          # Use -R or +R in 3j?
        EPRcons['negRlabel'] = False    # Use -R or +R in QN labels? (Will affect Xarray-based calculations.)
        EPRcons['negRcoordSwap'] = True    # Swap -R and +R in QN coords? (Will affect Xarray-based calculations.)

    elif phaseConvention == 'R':
        EPRcons['Rphase'] = True       # Apply (-1)^R phase?
        EPRcons['negR'] = True          # Use -R or +R in 3j?
        EPRcons['negRlabel'] = False    # Use -R or +R in QN labels? (Will affect Xarray-based calculations.)
        EPRcons['negRcoordSwap'] = True    # Swap -R and +R in QN coords? (Will affect Xarray-based calculations.)

    elif phaseConvention == 'E':
        EPRcons['Rphase'] = True        # Apply (-1)^R phase?
        EPRcons['negR'] = True          # Use -R or +R in 3j?
        EPRcons['negRlabel'] = False    # Use -R or +R in QN labels? (Will affect Xarray-based calculations.)
        EPRcons['negRcoordSwap'] = True    # Swap -R and +R in QN coords? (Will affect Xarray-based calculations.)

    phaseCons['EPR'] = EPRcons

    #*** For Lambda term (as set by MFproj())
    lambdaCons = {}
    if phaseConvention == 'S':
        lambdaCons['negMup'] = True     # Use -mup term in 3j?
        lambdaCons['negRp'] = True     # Use -Rp term in 3j?
        lambdaCons['phaseNegR'] = True  # Set for (-Rp, -R) phase convention (in Wigner D term)?
        lambdaCons['conjFlag'] = False  # Set for conjuate Wigner D?
        lambdaCons['RpPhase'] = True    # Apply (-1)^Rp phase term?

    elif phaseConvention == 'R':
        lambdaCons['negMup'] = True     # Use -mup term in 3j?
        lambdaCons['negRp'] = True     # Use -Rp term in 3j?
        lambdaCons['phaseNegR'] = True  # Set for (-Rp, -R) phase convention (in Wigner D term)?
        lambdaCons['conjFlag'] = False  # Set for conjuate Wigner D?
        lambdaCons['RpPhase'] = True    # Apply (-1)^Rp phase term?

    elif phaseConvention == 'E':
        lambdaCons['negMup'] = True     # Use -mup term in 3j?
        lambdaCons['negRp'] = True     # Use -Rp term in 3j?
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
        betaCons['negM'] = genMatEcons['negm']       # Use -M term in 3j? Should be anti-correlated with genMatEcons['negm']...? 31/03/20 NOW correlated with mfblmCons['Mphase']
        # betaCons['negM'] = not(genMatEcons['negm'])     # Use -M term in 3j? Should be anti-correlated with genMatEcons['negm']...? 31/03/20 NOW correlated with mfblmCons['Mphase']  18/07/23 - now testing after afblm phase tests below. THIS KILLS M!=0 for BLM term.
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

        mfblmCons['mupPhase'] = False            # Apply (-1)^(mup - p) phase term. Already incoporated into MFproj()?

        mfblmCons['BLMmPhase'] = False          # TESTING ONLY - switch signs (m, M) terms before 3j calcs.

    if phaseConvention == 'E':
        mfblmCons['negRcoordSwap'] = False       # Swap -R and +R in EPR Xarray coords? Already in EPRCons()?

        mfblmCons['negMcoordSwap'] = False       # Swap +/-M coords.
        mfblmCons['Mphase'] = betaCons['negM']   # Apply (-1)^M phase term. Correlated with +/-M term switch.

        mfblmCons['negmCoordSwap'] = False       # Swap +/-m coords.
        mfblmCons['mPhase'] = False              # Apply (-1)^m phase term.

        mfblmCons['mupPhase'] = True            # Apply (-1)^(mup - p) phase term. Already incoporated into MFproj()?

        mfblmCons['BLMmPhase'] = False          # TESTING ONLY - switch signs (m, M) terms before 3j calcs.

    phaseCons['mfblmCons'] = mfblmCons


    #*** For AFPAD product case, as calculated in afblmXprod()
    phaseCons['afblmCons'] = {}

    if phaseConvention == 'S':
        # (+/-)M phase selection, set as per existing code, betaCons['negM'] = genMatEcons['negm']       # Use -M term in 3j? Should be anti-correlated with genMatEcons['negm']...? 31/03/20 NOW correlated with mfblmCons['Mphase']
        # Note this is correlated with QN generation in genllpMatE() - should set equivalent fn for alignment terms.
        # In existing case this arises from M = (-m+mp) or M = -(m+mp) choice.
        phaseCons['afblmCons']['negM'] = phaseCons['genMatEcons']['negm']  # 01/07/23 - set to not() here, should be correct and tested for 3D alignment case (doesn't affect 1D results)
                                                                           # 17/07/23 - Added separate 'S' case, since for normal case this removes M!=0 terms from DeltaTerm!
                                                                           #            THIS MAY ALSO BE AN ISSUE FOR 'E' below, may have inconsistencies elsewhere?
        phaseCons['afblmCons']['negQ'] = True
        phaseCons['afblmCons']['negS'] = True

        phaseCons['afblmCons']['llpPhase'] = True  # Apply (-1)^(l-lp) phase term?

    if phaseConvention == 'R':
        # (+/-)M phase selection, set as per existing code, betaCons['negM'] = genMatEcons['negm']       # Use -M term in 3j? Should be anti-correlated with genMatEcons['negm']...? 31/03/20 NOW correlated with mfblmCons['Mphase']
        # Note this is correlated with QN generation in genllpMatE() - should set equivalent fn for alignment terms.
        # In existing case this arises from M = (-m+mp) or M = -(m+mp) choice.
        phaseCons['afblmCons']['negM'] = phaseCons['genMatEcons']['negm']  # 01/07/23 - set to not() here, should be correct and tested for 3D alignment case (doesn't affect 1D results)
                                                                           # 17/07/23 - Added separate 'S' case, since for normal case this removes M!=0 terms from DeltaTerm!
                                                                           #            THIS MAY ALSO BE AN ISSUE FOR 'E' below, may have inconsistencies elsewhere?
        phaseCons['afblmCons']['negQ'] = True
        phaseCons['afblmCons']['negS'] = True

        phaseCons['afblmCons']['llpPhase'] = True  # Apply (-1)^(l-lp) phase term?


    if phaseConvention == 'E':
        # (+/-)M phase selection, set as per existing code, betaCons['negM'] = genMatEcons['negm']       # Use -M term in 3j? Should be anti-correlated with genMatEcons['negm']...? 31/03/20 NOW correlated with mfblmCons['Mphase']
        # Note this is correlated with QN generation in genllpMatE() - should set equivalent fn for alignment terms.
        # In existing case this arises from M = (-m+mp) or M = -(m+mp) choice.
        phaseCons['afblmCons']['negM'] = not(phaseCons['genMatEcons']['negm'])  # 01/07/23 - set to not() here, should be correct and tested for 3D alignment case (doesn't affect 1D results)
        phaseCons['afblmCons']['negQ'] = True
        phaseCons['afblmCons']['negS'] = True

        phaseCons['afblmCons']['llpPhase'] = True  # Apply (-1)^(l-lp) phase term?


    #*** For LF testing with CG terms.
    phaseCons['lfblmCGCons'] = {}

    # l,m terms
    phaseCons['lfblmCGCons']['negmp'] = False # mp > -mp  sign flip
                    # Including this restricts things to mp=0 only? Phase con choice? Doesn't seem correct!
                    # Full calculation including this term sends PU continuum to zero/NaN.
    phaseCons['lfblmCGCons']['negM'] = True  # M > -M sign flip.

    # Photon terms
    phaseCons['lfblmCGCons']['negmup'] = False # mup > -mup sign flip.  See note above - kills many terms if True (incorrectly it seems)
    phaseCons['lfblmCGCons']['negMP'] = False  # M > -M sign flip, photon term.
                                            # Originally thought this was a typo, but do need to set False as per manuscript!


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
        - 'sel' : use full selection routine, :py:func:`epsproc.selQNsRow()`, should be most robust.
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

    # Int conversion UNLESS half int values!
    if not (uniquellpL%1.0).any():
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




#*************************************************************
#************** Core functions (Wigner/CG calculators)
#*************************************************************

# Tabulate Wigner 3j terms for a given problem/set of QNs
def w3jTable(Lmin = 0, Lmax = 10, QNs = None, mFlag = True, nonzeroFlag = False, form = '2d', dlist = ['l','lp','L','m','mp','M'], backend = 'par', verbose = 0):
    r"""
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
        NOTE: some return types will convert QNs to int when constructing Xarray output, unless half-int values present.
        Functions using :py:func:`epsproc.geomFunc.geomCalc.remapllpL()` support half-int values in Xarray.

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
    if backend == 'par':
        w3j_QNs = w3jprange(QNs)

    elif backend == 'sympy':
        w3j_QNs = w3jSympy(QNs)
        w3j_QNs = np.array(w3j_QNs)  # Convert to np.array to match other function returns

    else:
        w3j_QNs = Wigner3jQNs(QNs)  # Fallback to vanilla spherical_functions version.

    if nonzeroFlag:
        # Drop zero terms
        nzInd = np.nonzero(w3j_QNs)
        w3j_QNs = w3j_QNs[nzInd]
        QNs = QNs[nzInd]

    if form == '2d':

        return np.c_[QNs, w3j_QNs]


    elif form == 'xarray':
        # Multindex case - this retains same size as original, since coords are stacked

        # Int conversion UNLESS half int values!
        if not (QNs.real%1.0).any():
            QNs = pd.MultiIndex.from_arrays(QNs.real.T.astype('int8'), names = dlist)
        else:
            QNs = pd.MultiIndex.from_arrays(QNs.real.T, names = dlist)

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
        # QNs = pd.MultiIndex.from_arrays(QNs.real.T.astype('int8'), names = dlist)
        # Int conversion UNLESS half int values!
        if not (QNs.real%1.0).any():
            QNs = pd.MultiIndex.from_arrays(QNs.real.T.astype('int8'), names = dlist)
        else:
            QNs = pd.MultiIndex.from_arrays(QNs.real.T, names = dlist)

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



# Set CG function
def CG(QNs, dlist = ['l', 'lp', 'L', 'm', 'mp', 'M'], form = 'xarray'):
    """
    Basic Clebsch-Gordan from 3j calculation, from table of input QNs (corresponding to CG term defn.).

    This implements numerical defn. from Moble's Spherical Functions, https://github.com/moble/spherical_functions/blob/master/spherical_functions/recursions/wigner3j.py

    `def clebsch_gordan(j_1, m_1, j_2, m_2, j_3, m_3)`

    `(-1.)**(j_1-j_2+m_3) * math.sqrt(2*j_3+1) * Wigner3j(j_1, j_2, j_3, m_1, m_2, -m_3)`

    22/06/20 - barebones version for quick testing, should upgrade as per w3jTable (which is the back-end here in any case).

    Parameters
    ----------
    QNs : np.array
        List of QNs [l, lp, L, m, mp, M] to compute 3j terms for.

    form : string, optional, default = 'xarray'
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

    Returns
    -------


    """

    # Set phase convention, CG(M3) > 3j(-M3)
    QNs[:,5] *= -1

    # Generate some 3j values (Xarray form) from supplied QNs
    w3j = w3jTable(QNs = QNs, dlist = dlist, form = form, nonzeroFlag = True)

    # Phase and degen terms
    # For testing set explicit version for full (l,lp,L) term, and photon (1,1,L) term
    # if 'l' in dlist:
    #     CGphase = np.power(-1, np.abs(w3j.l - w3j.lp + w3j.M))
    #     # CGphase = np.power(-1, np.abs(w3j.lp - w3j.l + w3j.M))  # With alternative ordering (lp,l)
    # else:
    #     CGphase = np.power(-1, np.abs(w3j.M))
    #
    # CGdegen = np.sqrt(2*w3j.L + 1)

    # Replaces above - use dim labels as passed!
    CGphase = np.power(-1, np.abs(w3j[dlist[0]] - w3j[dlist[1]] + w3j[dlist[5]]))
    CGdegen = np.sqrt(2*w3j[dlist[2]] + 1)


    return CGphase * CGdegen * w3j



#*************************************************************
#************** Geometric terms (EPR, betaTerm etc.)
#*************************************************************

def EPR(QNs = None, p = None, ep = None, normE = False,
        nonzeroFlag = True, form = '2d', dlist = None,
        phaseConvention = 'S', verbose = 0):
    r"""Define polarization tensor (LF) for 1-photon case.

    Define field terms (from QM book, corrected vs. original S\&U version
    - see ``beta-general-forms\_rewrite\_290917.lyx``):

    .. math::
        \begin{equation}
        E_{PR}(\hat{e})=[e\otimes e^{*}]_{R}^{P}=[P]^{\frac{1}{2}}\sum_{p}(-1)^{R}\left(\begin{array}{ccc}
        1 & 1 & P\\
        p & R-p & -R
        \end{array}\right)e_{p}e_{R-p}^{*}
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

    normE : bool, optional, default = None
        If True, renorm ep field strength terms to unity.

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

    verbose : optional, default=0
        Set verbosity.

    Examples
    ---------
    # Generate full EPR list with defaults
    >>> EPRtable = EPR()

    # Return as Xarray
    >>> EPRtable = EPR(form = 'xarray')

    Note
    ----

    18/07/23: fixed issues with dlist > Xarray conversion, and switched to setting in function to avoid odd behaviour with passed list.

    01/03/24: implemented ep terms for 2d and xarray forms.
    Default case with ep = np.ones(len(p)) should match previous version outputs (==all terms).
    TODO: test with new pol functions (see ep.efields.epol).

    """
    # 18/07/23: use dlist=None as default.
    # This fixes issues with repeated calls, see https://stackoverflow.com/questions/1132941/least-astonishment-and-the-mutable-default-argument
    if dlist is None:
        dlist = ['l', 'lp', 'P', 'p', 'R-p', 'R']

    if verbose >1:
        print(locals())

    # Set phase conventions
    phaseCons = setPhaseConventions(phaseConvention = phaseConvention)

    if verbose >2:
        print(phaseCons)

    if phaseCons['EPR']['negRlabel']:
        dlist[-1] = '-R'    # Set -R in QN labels.  Note this assumes dlist[-1] = 'R'

    if verbose:
        print(f"Set dlist={dlist}")

    # Set dim labels for reference/Xarray case - now passed
    # dlist = ['l', 'lp', 'P', 'p', 'R-p', 'R']

    # If p is not passed, set for all allowed values
    if p is None:
        p = [-1,0,1]

    # If no field strength components are defined, set to unity
    if ep is None:
        ep = np.ones(len(p))

    # Force to numpy, otherwise indexing later may fail
    if type(ep) is not np.ndarray:
        ep = np.array(ep)

    if normE:
        # Norm field?
        ep = ep/np.sqrt((ep**2).sum())

    if form == 'xarray':
        epXR = xr.DataArray(ep, coords={'p':p}, dims='p')
        # Conj term, p > R-p for cross-product
        eRpXR = epXR.conj().copy().rename({'p':'R-p'})

        # Include additional phase term?
        # TODO: consider this, see Zare p209.
        # Only required for case with real fields defined, and cos+sin components.
        # RpPhase = (-1)**eRpXR.coords['R-p'].pipe(np.abs)

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
    EPRtable = None   # Testing 18/07/23, set None explictly here, seeing some weird caching issues when testing.
    EPRtable = w3jTable(QNs = QNs, nonzeroFlag = nonzeroFlag, form = form, dlist = dlist)
    # EPRtable = w3jTable(QNs = QNs, nonzeroFlag = nonzeroFlag, form = form, dlist = dlist, backend='default')  # Testing 18/07/23, set default backend here, seeing some weird caching issues when testing.

    # Phase & degen terms
    if form == '2d':
        Rphase = np.power(-1, np.abs(EPRtable[:,5]))
        Pdegen = np.sqrt(2*EPRtable[:,2] + 1)
        EPRtable[:,-1] = EPRtable[:,-1]*Rphase*Pdegen

        # Set Ep field strength terms and multiply
        pRPInd = EPRtable[:,[3,4]]  # Select (p, R-p) terms only

        # Multiplication vectors from indexes into ep values...
        epMult = ep[(pRPInd[:,0]+1).astype(int)]
        epMultConj = ep[(pRPInd[:,1]+1).astype(int)].conj()

#         Test terms - resultant is 0 if ep=p
#         pRPInd[:,0] - epMult
#         pRPInd[:,1] - epMultConj

        # Update main table
        EPRtable[:,-1] = EPRtable[:,-1] * epMult * epMultConj


    elif form == 'xarray':
        Pdegen = np.sqrt(2*EPRtable.P + 1)  # Note sqrt here - as per U&S defn.
        # Pdegen = 2*EPRtable.P + 1

        # Set phase choice for 3j term
        if phaseCons['EPR']['Rphase']:
            # Rphase = np.power(-1, np.abs(EPRtable.R))
            Rphase = np.power(-1, np.abs(EPRtable[dlist[-1]]))  # 18/07/23 - fixed for var name
            EPRtable *= Rphase*Pdegen
            # EPRtable = EPRtable*Rphase*Pdegen
        else:
            EPRtable *= Pdegen
            # EPRtable = EPRtable*Pdegen

        # Switch coord sign?
        if phaseCons['EPR']['negRcoordSwap']:
            temp = EPRtable.unstack()   # Need to unstack to change MultiIndex coords?  Might be a cleaner way?
            # temp[dlist[-1]] *= -1       # Label from dlist to allow for +/-R label here.
            temp[dlist[-1]] = -1 * temp[dlist[-1]]       # 18/07/23 - updated for XR > v0.15
            EPRtable = temp.stack({'QN':dlist}).dropna(dim = 'QN',how = 'all')
            # NOTE: Without dropna here dims grow! Default settings have 18 elements, but end up with 135 and lots of NaNs.

        # Multiply by field strength array
        # TO TEST: does QN coord swap propagate correctly here?
        # Should be OK since ep terms just indexed by (p,p*)?
        EPRmult = EPRtable.unstack() * epXR * eRpXR
        EPRtable = EPRmult.stack({'QN':dlist}).dropna(dim = 'QN',how = 'all')

        EPRtable.attrs['dataType'] = 'EPR'
        EPRtable.attrs['phaseCons'] = phaseCons

    return EPRtable


# BetaTerm for BLM tensor coupling.
def betaTerm(QNs = None, Lmin = 0, Lmax = 10, nonzeroFlag = True, form = '2d', dlist = ['l', 'lp', 'L', 'm', 'mp', 'M'], phaseConvention = 'S'):
    r"""Define BLM coupling tensor

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

    # Set QNs, if not supplied - run here (rather than in w3jTable) to provide local copy
    if QNs is None:
        QNs = genllL(Lmin = Lmin, Lmax = Lmax, mFlag = True)

    # Set -M for 3j term if required
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
        # try:
        #     # 3j product term
        #     # BLMtable *= mPhase*np.sqrt(degen)*BLMtable.sel(m=0,mp=0,M=0)  # Accidental correlations in results...?
        #     BLMtable = BLMtable*BLMtable.sel(m=0,mp=0,M=0).drop('mSet').squeeze()  # NOTE - drop dims here to prevent correlated product.
        #
        # # If (0,0,0) terms are not already calculated, do so.
        # # Might be cleaner just to use this in all cases?
        # except KeyError:
        #     # Case for using matE directly, set for m=0 terms only.
        #     # thrj0 = w3jTable(QNs = genllpMatE(matE, mFlag = False), nonzeroFlag = True, form = form, dlist = dlist)
        #     # Case for list of QNs
        #     thrj0 = w3jTable(QNs = genllLList(QNs, uniqueFlag = True, mFlag = False), nonzeroFlag = True, form = form, dlist = dlist)
        #     BLMtable = BLMtable*thrj0.drop('mSet').squeeze()

        # This version currently doesn't work - WHY?  Seems to drop necessary (l,lp,L) correlations?
        # thrj0 = w3jTable(QNs = genllLList(QNs, uniqueFlag = True, mFlag = False), nonzeroFlag = True, form = form, dlist = dlist)
        # BLMtable = BLMtable*thrj0.drop('mSet').squeeze()

        # Now with unstack to ensure QN dims.
        # TODO: this now breaks xdaLM expected output (stacked) - needs rationalisation over all fns however!
        # Pass also local QNs to ensure phase conventions set above.
        thrj0 = w3jTable(QNs = genllLList(QNs, uniqueFlag = True, mFlag = False), nonzeroFlag = True, form = form, dlist = dlist)
        BLMtable = BLMtable.unstack()*thrj0.drop('mSet').squeeze().unstack()

        mPhase = np.power(-1, np.abs(BLMtable.m))
        degen = (2*BLMtable.l+1)*(2*BLMtable.lp+1)*((2*BLMtable.L+1))/(4*np.pi)

        if phaseCons['betaCons']['mPhase']:
            BLMtable *= mPhase*np.sqrt(degen)
            # BLMtable = BLMtable*mPhase*np.sqrt(degen)
        else:
            BLMtable *= np.sqrt(degen)
            # BLMtable = BLMtable*np.sqrt(degen)

        BLMtable.attrs['dataType'] = 'betaTerm'
        BLMtable.attrs['phaseCons'] = phaseCons

    else:
        print(f"Form {form} not implemented.")
        BLMtable = None

    return BLMtable


# Define lambdaTerm, MF projection.
def MFproj(QNs = None, RX = None, nonzeroFlag = True, form = '2d', dlist = ['l', 'lp', 'P', 'mu', 'mup', 'Rp', 'R'],
            eNames = ['Ph','Th','Ch'], phaseConvention = 'S'):
    r"""
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

    eNames : optional, list, default = ['Ph','Th','Ch']
        Set names for Euler angles in output.
        Note:

        - eNames = ['P','T','C'] matches defaults in :py:func:`epsproc.sphCalcs.setPolGeoms()` and :py:func:`epsproc.sphCalcs.wDcalc()`, but conflicts with QNs 'P'.
        - eNames = ['Phi','Theta','Chi'] may give issues elsewhere, e.g. when multiplying by Ylms with same dim names.


    Notes
    -----
    This is very similar to :math:`E_{PR}` term.

    12/08/22    Added option for Euler angle dim names to avoid conflicts elsewhere.

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
                    # for Rp in np.arange(-P, P+1):  # Allow all Rp
                    # Rp = -(mu+mup)   # Fix Rp terms - not valid here, depends on other phase cons!
                    # for R in np.arange(-P, P+1):
                    #     # QNs.append([l, lp, P, mu, -mup, R, Rp])
                    #     if phaseCons['lambdaCons']['negRp']:
                    #         Rp *= -1
                    #     if phaseCons['lambdaCons']['negMup']:
                    #         QNs.append([l, lp, P, mu, -mup, Rp, R])   # 31/03/20: FIXED bug, (R,Rp) previously misordered!!!
                    #     else:
                    #         QNs.append([l, lp, P, mu, mup, Rp, R])

                    # Rearranged for specified Rp case
                    for R in np.arange(-P, P+1):
                        # if phaseCons['lambdaCons']['negMup']:
                        #     mup = -mup

                        if phaseCons['lambdaCons']['negRp']:
                            # Rp = mu+mup
                            Rp = mup - mu
                        else:
                            Rp = -(mu+mup)

                        # Switch mup sign for 3j?  To match old numerics, this is *after* Rp assignment (sigh).
                        if phaseCons['lambdaCons']['negMup']:
                            mup = -mup

                        QNs.append([l, lp, P, mu, mup, Rp, R])

        QNs = np.array(QNs)

    # Set to calculate for default (x,y,z) pol geoms.
    if RX is None:
        RX = setPolGeoms()

    #*********************** Calculate terms
    # 3j term with w3jTable()
    # lambdaTable = w3jTable(QNs = QNs[:,:-1], nonzeroFlag = nonzeroFlag, form = form, dlist = dlist[:-1])  # Pass only dims for 3j term
    if form == '2d':
        nonzeroFlag = False   # For 2d case override to ensure consistent size.

    lambdaTable = w3jTable(QNs = QNs, nonzeroFlag = nonzeroFlag, form = form, dlist = dlist)
    # lambdaTable = w3jTable(QNs = QNs[:,0:5], nonzeroFlag = nonzeroFlag, form = form, dlist = dlist[0:5])  # Pass all QNs to keep R label. w3jpRange will select/use QN[:,0:5] only... but may have issues with Xarray dims!!!

    # D term with wDcalc(), includes order switch for (R,Rp)
    # Pass RX.data since wDcalc is curently only set for np.array inputs.
    # Set XFlag by form.
    if form.startswith('x'):
        XFlag = True
    else:
        XFlag = False

    # Subselect on QNs and calculate
    QNind = np.array([2, 5, 6])   # Set (P,Rp,R)
    dRed = [dlist[n] for n in QNind]

    QNwD = QNs[:,QNind]
    QNun = np.unique(QNs[:,QNind], axis=0)

    # Set for (-Rp, -R) phase convention
    if phaseCons['lambdaCons']['phaseNegR']:
        QNwD[:,1:] = -1 * QNwD[:,1:]
        QNun[:,1:] = -1 * QNun[:,1:]

    # Correct for case where -Rp is already set in 3j
    # ONLY want to do this is setting phase convention explicitly here?
    if phaseCons['lambdaCons']['negRp']:
        QNwD[:,1] = -1 * QNwD[:,1]
        QNun[:,1] = -1 * QNun[:,1]

    if form == '2d':
        # Cal for all values - in cases with duplicate QNs this may lead to indexing issues later in Xarray case.
        # Should be OK for 2d case however, will provide table matching full QN list.
        # NOTE in 2d case, wDcalc() also outputs R, QNs - skip these since they're already set in this case
        lambdaD, *RQN = wDcalc(QNs = QNwD, R = RX.data, XFlag = XFlag, dlist = dRed, eNames = eNames,
                                conjFlag = phaseCons['lambdaCons']['conjFlag'])
        lambdaD = np.asarray(lambdaD)

    elif form.startswith('x'):
        # Calc for unique values only to avoid duplicate coords

        # Special case for single pol geom, otherwise passes 0D array with a single value and causes issues later.
        if RX.size == 1:
            lambdaD = wDcalc(QNs = QNun, R = RX.data.item(), XFlag = XFlag, dlist = dRed, eNames = eNames,
                             conjFlag = phaseCons['lambdaCons']['conjFlag'])

            lambdaD['Labels']=('Euler', [RX.Labels.item()])  # Propagate labels, currently wDcalc only takes RX.data

        else:
            lambdaD = wDcalc(QNs = QNun, R = RX.data, XFlag = XFlag, dlist = dRed, eNames = eNames,
                             conjFlag = phaseCons['lambdaCons']['conjFlag'])

            # lambdaD['Labels']=('Euler', RX.Labels.values)  # Propagate labels, currently wDcalc only takes RX.data
            # lambdaD['Labels']=('Euler', RX.Labels)  # Propagate labels, currently wDcalc only takes RX.data
            lambdaD['Labels']=('Euler', RX.Labels.data)  # 09/06/22 reinstated this during testing, in XR v2022.3.0 the above line gives "TypeError: Using a DataArray object to construct a variable is ambiguous, please extract the data using the .data property."
                                                         # Should be back-compatible? Or change to .values as previously?

        lambdaD = lambdaD.swap_dims({'Euler':'Labels'})  # Swap dims to labels.
        lambdaD.attrs['dataType'] = 'WignerD'



    #***************** Multiplications & phase
    if form.startswith('x'):
        lun = lambdaTable.unstack('QN')
        lDun = lambdaD.unstack('QN')

        # Reset phase choices here to allow for correct multiplication of +Rp and D(P,-Rp,-R) terms in Xarray
        if phaseCons['lambdaCons']['phaseNegR']:
            lDun['R'] = -1 * lDun['R']
            lDun['Rp'] = -1 * lDun['Rp']

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
            Rpphase = np.power(-1, np.abs(QNs[:,-2]))  # 08/04/20 fixed indexing for Rp (not R!)
        else:
            Rpphase = np.ones(QNs[:,-2].size)

        # Loop over sets of Euler angles
        lTerm = []
        for eInd in range(lambdaD.shape[1]):
            lTerm.append(Rpphase * lambdaD[:,eInd] * lambdaTable[:,-1])

        lTerm = np.array(lTerm).T

    else:
        print(f'Form {form} not supported.')


    return lTerm, lambdaTable, lambdaD, QNs



# Calculate alignment term - this cell should form core function, cf. betaTerm() etc.
def deltaLMKQS(EPRX, AKQS, phaseConvention = 'S'):
    r"""
    Calculate aligned-frame "alignment" term:

    .. math::
        \begin{equation}
        \sum_{K,Q,S}\Delta_{L,M}(K,Q,S)A_{Q,S}^{K}(t)
        \end{equation}

    .. math::
        \begin{equation}
        \Delta_{L,M}(K,Q,S)=(2K+1)^{1/2}(-1)^{K+Q}\left(\begin{array}{ccc}
        P & K & L\\
        R & -Q & -M
        \end{array}\right)\left(\begin{array}{ccc}
        P & K & L\\
        R' & -S & S-R'
        \end{array}\right)
        \end{equation}

    15/06/20 IN PROGRESS

    Parameters
    ----------
    EPRX : Xarray
        Polarization terms in an Xarray, as set by :py:func:`epsproc.geomCalc.EPR`

    AKQS : Xarray
        Alignement terms in an Xarray, as set by :py:func:`epsproc.setADMs`

    Returns
    -------
    AFterm : Xarray
        Full term, including multiplication and sum over (K,Q,S) (note S-Rp term is retained).

    DeltaKQS : Xarray
        Alignment term :math:`\Delta_{L,M}(K,Q,S)`.

    To do
    -----
    - Add optional inputs.
    - Add error checks.
    See other similar functions for schemes.

    """

    # Get phase conventions
    phaseCons = setPhaseConventions(phaseConvention = phaseConvention)

    # Set QNs
    QNs1, QNs2 = genKQStermsFromTensors(EPRX, AKQS, uniqueFlag = True, phaseConvention = phaseCons)

    # Then calc 3js.... as per betaTerm
    form = 'xdaLM'  # xds
    dlist1 = ['P', 'K', 'L', 'R', 'Q', 'M']
    dlist2 = ['P', 'K', 'L', 'Rp', 'S', 'S-Rp']

    # Copy QNs and apply any additional phase conventions
    QNs1DeltaTable = QNs1.copy()
    QNs2DeltaTable = QNs2.copy()

    # Set additional phase cons here - these will be set in master function eventually!
    # NOW - set in setPhaseConventions()
    # # NOTE - only testing for Q=S=0 case initially.
    # phaseCons['afblmCons']['negM'] = phaseCons['genMatEcons']['negm']  # IF SET TO TRUE THIS KNOCKS OUT M!=0 terms - not sure if this is correct here, depends also on phase cons in genKQStermsFromTensors().
    #                                                                     # Yeah, looks like phase error in current case, get terms with R=M, instead of R=-M
    #                                                                     # Confusion is due to explicit assignment of +/-M terms in QN generation (only allowed terms), which *already* enforces this phase convention.
    # phaseCons['afblmCons']['negQ'] = True
    # phaseCons['afblmCons']['negS'] = True

    # Switch sign Q > -Q before 3j calcs.
    if phaseCons['afblmCons']['negQ']:
        QNs1DeltaTable[:,4] *= -1

    # Switch signs M > -M before 3j calcs.
    if phaseCons['afblmCons']['negM']:
        QNs1DeltaTable[:,5] *= -1

    # Switch sign S > -S before 3j calcs.
    if phaseCons['afblmCons']['negS']:
        QNs2DeltaTable[:,4] *= -1


    # Calculate two 3j terms, with respective QN sets
    thrj1 = w3jTable(QNs = QNs1DeltaTable, nonzeroFlag = True, form = form, dlist = dlist1)
    thrj2 = w3jTable(QNs = QNs2DeltaTable, nonzeroFlag = True, form = form, dlist = dlist2)

    # Multiply
    thrjMult = thrj1.unstack() * thrj2.unstack()

    # Additional terms & multiplications
    Kdegen = np.sqrt(2*thrjMult.K + 1)
    KQphase = np.power(-1, np.abs(thrjMult.K + thrjMult.Q))

    DeltaKQS =  Kdegen * KQphase * thrjMult

    # AF term
    AFterm = (DeltaKQS * AKQS.unstack()).sum({'K','Q','S'})

    return AFterm, DeltaKQS

# Calculate alignment term - this cell should form core function, cf. betaTerm() etc.
# 13/01/21 Rough version for AF wavefunction expansion, adapted from existing deltaLMKQS() function.
def deltalmKQSwf(matE, AKQS, phaseConvention = 'S', dlist1 = ['l','E','K','m','mu','Q'], dlist2 = ['l','E','K','Lambda','mu0','S']):
    r"""
    Calculate aligned-frame "alignment" term, version for AF wavefunction expansion.

    .. math::
        \begin{equation}
        ^{AF}\Delta_{l,m}(K,Q,S)=(-1)^{m-\Lambda}(-1)^{\mu-\mu_{0}}(-1)^{Q-S}\left(\begin{array}{ccc}
        l & 1 & K\\
        -m & -\mu & -Q
        \end{array}\right)\left(\begin{array}{ccc}
        l & 1 & K\\
        -\Lambda & -\mu_{0} & -S
        \end{array}\right)
        \end{equation}

    13/01/21 IN PROGRESS, adapted from existing deltaLMKQS() function.

    NOTE: photon QN currently labelled as (E,mu,mu0), but may want to change to avoid confusion with EPR term.

    Parameters
    ----------
    matE : Xarray
        Xarray containing matrix elements, with QNs (l,m), as created by :py:func:`readMatEle`

    AKQS : Xarray
        Alignement terms in an Xarray, as set by :py:func:`epsproc.setADMs`

    phaseConvention : optional, str, default = 'S'
        Set phase conventions with :py:func:`epsproc.geomCalc.setPhaseConventions`.
        To use preset phase conventions, pass existing dictionary.
        If matE.attrs['phaseCons'] is already set, this will be used instead of passed args.

    dlist1, dlist2 : optional, lists, defaults =  ['l','E','K','m','mu','Q'], ['l','E','K','Lambda','mu0','S']
        Labels for output QNs.

    Returns
    -------
    AFterm : Xarray
        Full term, including multiplication and sum over (K,Q,S).

    DeltaKQS : Xarray
        Alignment term :math:`\Delta_{L,M}(K,Q,S)`.

    To do
    -----
    - Add optional inputs.
    - Add error checks.
    See other similar functions for schemes.

    """

    # Get phase conventions
    phaseCons = setPhaseConventions(phaseConvention = phaseConvention)

    # Set QNs
    QNs1, QNs2 = genLmLamKQStermsFromTensors(matE, AKQS, uniqueFlag = True, phaseConvention = phaseConvention)

    # Then calc 3js.... as per betaTerm
    form = 'xdaLM'  # xds
    # dlist1 = ['P', 'K', 'L', 'R', 'Q', 'M']  # Now passed to fn.
    # dlist2 = ['P', 'K', 'L', 'Rp', 'S', 'S-Rp']

    # Copy QNs and apply any additional phase conventions
    QNs1DeltaTable = QNs1.copy()
    QNs2DeltaTable = QNs2.copy()

    # Set additional phase cons here - these will be set in master function eventually!
    # NOW - set in setPhaseConventions()
    # # NOTE - only testing for Q=S=0 case initially.
    # phaseCons['afblmCons']['negM'] = phaseCons['genMatEcons']['negm']  # IF SET TO TRUE THIS KNOCKS OUT M!=0 terms - not sure if this is correct here, depends also on phase cons in genKQStermsFromTensors().
    #                                                                     # Yeah, looks like phase error in current case, get terms with R=M, instead of R=-M
    #                                                                     # Confusion is due to explicit assignment of +/-M terms in QN generation (only allowed terms), which *already* enforces this phase convention.
    # phaseCons['afblmCons']['negQ'] = True
    # phaseCons['afblmCons']['negS'] = True

    # Switch signs (m,M) before 3j calcs.
    # if phaseCons['afblmCons']['negQ']:
    #     QNs1DeltaTable[:,4] *= -1
    #
    # # Switch sign Q > -Q before 3j calcs.
    # if phaseCons['afblmCons']['negM']:
    #     QNs1DeltaTable[:,5] *= -1
    #
    # # Switch sign S > -S before 3j calcs.
    # if phaseCons['afblmCons']['negS']:
    #     QNs2DeltaTable[:,4] *= -1


    # Calculate two 3j terms, with respective QN sets
    thrj1 = w3jTable(QNs = QNs1DeltaTable, nonzeroFlag = True, form = form, dlist = dlist1)
    thrj2 = w3jTable(QNs = QNs2DeltaTable, nonzeroFlag = True, form = form, dlist = dlist2)

    # Multiply
    thrjMult = thrj1.unstack() * thrj2.unstack()

    # Additional terms & multiplications
    # NOTE THESE ASSUME QN labels.
    # Kdegen = np.sqrt(2*thrjMult.K + 1)
    # degen = 8*(np.pi**2)
    degen = 1
    lLamPhase = np.power(-1, np.abs(thrjMult.m - thrjMult.Lambda))
    muPhase = np.power(-1, np.abs(thrjMult.mu - thrjMult.mu0))
    QSphase = np.power(-1, np.abs(thrjMult.Q - thrjMult.S))

    DeltaKQS =  degen * lLamPhase * muPhase * QSphase * thrjMult

    # AF term
    AFterm = (DeltaKQS * AKQS.unstack()).sum({'K','Q','S'})

    return AFterm, DeltaKQS
