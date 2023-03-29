"""
ePSproc routines for Spherical Harmonic function conversions

- Convert Real > Complex spherical harmonics with sphRealConvert()
- Calculate conjugates with sphConj()
- Convert to and from SHtools objects with SHcoeffsFromXR() and XRcoeffsFromSH()

v2  29/03/23    Added basic "cleanLMcoords" routine, may need some work.
v1  20/09/22

"""

import numpy as np
import xarray as xr
import copy    # For attrs deepcopy.

from epsproc.util.listFuncs import genLM, YLMtype, YLMdimList
from epsproc.util.conversion import multiDimXrToPD
from epsproc.util.misc import setDefaultArgs

try:
    import pyshtools as pysh
except ImportError as e:
    if e.msg != "No module named 'pyshtools'":
        raise
    print('* pyshtools not found, SHtools functions not available. If required, run "pip install pyshtools" or "conda install -c conda-forge pyshtools" to install.')


# Set SHtools object and init with set_coeffs
# x.set_coeffs(values, ls, ms)
# NOTE THIS CASTS TO REAL for real harmonic, so not sure how to convert dipoles here - will only apply to real part...?

# npTab = dfLong.to_numpy()
# clm = pysh.SHCoeffs.from_zeros(lmax = 6)   # Defaults to kind = 'real'

# def SHcoeffsFromXR(dataIn, kind = None, keyDims = None):
def SHcoeffsFromXR(dataIn, keyDims = None, **kwargs):
    """
    Xarray Spherical Harmonic coeffs to SHtools

    Parameters
    ----------
    dataIn : Xarray
        Data to convert.

    keyDims : dict, optional, default = None
        Passed to checkSphDims(dataIn, keyDims).
        If None use defaults.

    **kwargs : optional conversion args.
        Any valid args for conversion.
        If not specified dataIn.atts['harmonics'] settings will be used present, otherwise defaults will be used.
        Note names may not match SHtools naming.
        Defaults set as:
        defaults = {'kind':'complex','normType':'ortho','csPhase':True,}


    Notes
    ------

    MAY HAVE BETTER VERSION ELSEWHERE?

    This will only work for 1D array -- will need to slice/group for more general handling.

    TODO: update with general l,m dim names and handling, see checkSphDims()

    """

    # # Use settings on input array, or passed kind if set
    # if ('kind' in dataIn.attrs.keys()) and (kind is None):
    #     kind = dataIn.attrs['kind']
    #
    # elif ('harmonics' in dataIn.attrs.keys()) and (kind is None):
    #     kind = dataIn.attrs['harmonics']['kind']
    #
    # elif kind is None:
    #     kind = 'complex'

    # UPDATE 28/09/22 - use flexible setDefaultArgs function.
    defaults = {'kind':'complex',
                # 'keyDims':ep.listFuncs.YLMdimList(sType = 'sDict'),     # Set to default dim names.}
                'normType':'ortho',
                'csPhase':True,
               }

    if 'harmonics' in dataIn.attrs.keys():
        presetDict = dataIn.attrs['harmonics']

    setDefaultArgs(**locals(), **kwargs)  # Use setDefaults to check defaults vs. presetDict and any other locals()
                                          # NEED **kwargs to prevent nested dicts in locals()? Must be a better way?
                                          # Note this updated defaults dict in place.


    # Dim handling
    # TODO: may want to copy here? Should always be OK though.
    dataIn = checkSphDims(dataIn, keyDims)
    mDim = dataIn.attrs['harmonics']['mDim']  # For brevity below.
    lDim = dataIn.attrs['harmonics']['lDim']  # For brevity below.

#     else:
#         kind = kind

    # Init zeros SHtools coeffs object & populate
    # clm = pysh.SHCoeffs.from_zeros(lmax = dataIn[lDim].max().values, kind = kind)  # Minimal call

    # Include defaults{}. Note can't simply use **defaults as some names/items are different!
    clm = pysh.SHCoeffs.from_zeros(lmax = dataIn[lDim].max().values, kind = defaults['kind'], normalization = defaults['normType'], csphase = 1 if defaults['csPhase'] else -1)
    clm.set_coeffs(dataIn.values, dataIn[lDim].astype(int), dataIn[mDim].astype(int))  # NEEDS (values, ls, ms)

    return clm




def XRcoeffsFromSH(dataIn, keepSH = True, keyDims = None):
    """
    Convert SHtools coeffs object to Xarray

    NOTE: see ep.genLM and relate utils.
    ep.sphCalc.setBLMs  is similar too, although from list not SHtools object.

    Parameters
    -----------

    dataIn : SHtools object
        Data to be converted to Xarray.

    keepSH : bool, optional, default = True
        Set SHtools object to xr.attrs['SH'] if true.

    """

    # if keyDims is None:
    #     keyDims = YLMdimList(sType = 'sDict')   # Set to default dim names.

    # Set from inds
#     # Generate indexes
#     lmax = dataIn.lmax
#     LMlist = genLM(lmax).T
#     ind = pd.MultiIndex.from_arrays(LMlist, names = ['l','m'])

#     # Set +ve and -ve arrays to match SHtools convention [+/-,l,m] index
#     posMind = LMlist[1,:] >-1
#     negMind = LMlist[1,:] <1
#     indPos = pd.MultiIndex.from_arrays(LMlist[:,posMind], names = ['l','m'])
#     indNeg = pd.MultiIndex.from_arrays(LMlist[:,negMind], names = ['l','m'])


#     return indPos


#     # Direct to XR with sign
#     import numpy as np
#     lmXRTest = xr.DataArray(clmC.coeffs, coords={'sign':[0,1],'l':np.arange(0,clmC.coeffs.shape[1]), 'm':np.arange(0,clmC.coeffs.shape[2])})   #, coords={'l':})
#     lmXRTest


    # Split and restack - might be easier?
    lmXRTest1 = xr.DataArray(dataIn.coeffs[0,:,:], coords={'l':np.arange(0,dataIn.coeffs.shape[1]), 'm':np.arange(0,dataIn.coeffs.shape[2])})   #, coords={'l':})
    lmXRTest2 = xr.DataArray(dataIn.coeffs[1,:,:], coords={'l':np.arange(0,dataIn.coeffs.shape[1]), 'm': -1*np.arange(0,dataIn.coeffs.shape[2])})   #, coords={'l':})

    # Issues with duplicate l,m=0,0 case here
    # lmXRTest = xr.concat([lmXRTest1,lmXRTest2], dim='m', compat='no_conflicts', join='override')

    # Gives NaNs for duplicate case
    lmXRTest = xr.concat([lmXRTest1,lmXRTest2.where(lmXRTest2.l>0).where(lmXRTest2.m<0)], dim='m')  #, compat='no_conflicts', join='override')

    # Clean up - HAVE CODE FOR THIS ELSEWHERE?
    lmXRTestClean = lmXRTest.where(np.abs(lmXRTest.m)<=lmXRTest.l, drop=True)
    lmXRTestClean = lmXRTestClean.stack({'LM':['l','m']}).dropna(dim='LM')

    # Threshold?
    lmXRTestClean = lmXRTestClean.where(np.abs(lmXRTestClean)>1e-4).dropna(dim='LM')

    lmXRTestClean.attrs['SH'] = dataIn

    # Set harmonic attrs
    # lmXRTestClean.attrs['kind'] = dataIn.kind
    lmXRTestClean.attrs['harmonics'] = YLMtype(dtype = 'Harmonics from SHtools', kind = dataIn.kind, normType = dataIn.normalization,
                                                csPhase=dataIn.csphase, header = dataIn.header)

    # Set other attribs
    lmXRTestClean.attrs['harmonics'].update({'method': {'SHtools':'XRcoeffsFromSH'},
                                                'conj':None,
                                                'keyDims':keyDims,
                                                'Lrange':[0, dataIn.lmax],
                                                'res':None,
                                                'convention':None})

    return lmXRTestClean



def checkSphDims(dataIn, keyDims = None):
    """
    Very basic dim handling/checking for spherical harmonic Xarrays.

    Parameters
    ----------
    dataIn : Xarray
        Xarray containing values to convert.
        Must contain keyDims (stacked or unstacked).
        TODO: generalise this, should get from :py:func:`epsproc.util.listFuncs.getRefDims()`.

    keyDims : dict, optional, default = {'LM':['l','m']}
        Defines dim names for harmonic rank and order, and stacked dim name.


    Returns
    -------
    dataIn : Xarray
        As dataIn, but with additional .attrs['harmonics'] defining dim relations, and unstacked along keyDims[stackDim].


    TODO: should use ep.util.misc.checkDims here?
    TODO: check other dim handling functionality, may be reinventing things here.
    TODO: flag or separate func for unstack? May be confusing otherwise.
    TODO: check .attrs['harmonics']? This is currently implemented in unstackSphDims(), but not here.

    """

    if keyDims is None:
        keyDims = YLMdimList(sType = 'sDict')   # Set to default dim names.

        if 'BLM' in dataIn.dims:
            # dataIn = dataIn.rename({'BLM':'LM'})    # Switch naming to match defaults - UGLY.
            keyDims['BLM'] = keyDims.pop('LM')


    LMStackFlag = False
    stackDim = list(keyDims.keys())[0]
    dimList = keyDims[stackDim]
    lDim = keyDims[stackDim][0]
    mDim = keyDims[stackDim][1]

    if stackDim in dataIn.dims:
        # dataIn =  dataIn.unstack(stackDim)
        LMStackFlag = True

    if 'harmonics' not in dataIn.attrs.keys():
        dataIn.attrs['harmonics'] = {}

    dataIn.attrs['harmonics'].update({k:v for k,v in locals().items() if k != 'dataIn'})

    return dataIn


def unstackSphDims(dataIn, keyDims = None):
    """
    Check and unstack spherical harmonic (key) dims.

    TODO: handle stacking also?
    """

    if ('harmonics' not in dataIn.attrs.keys()) or ('stackDim' not in dataIn.attrs['harmonics'].keys()):
        dataIn = checkSphDims(dataIn, keyDims)

    if dataIn.attrs['harmonics']['stackDim'] in dataIn.dims:
        dataIn =  dataIn.unstack(dataIn.attrs['harmonics']['stackDim'])
        dataIn.attrs['harmonics']['LMStackFlag'] = True

    return dataIn


def tabulateLM(data, fillna = True):
    """
    Display set of B(L,M) values using Pandas.

    Thin wrapper for :py:func:`epsproc.multiDimXrToPD` for BLM data types or SHtools object.

    NOTE: currently returns pdTab, may want to switch to display.
    """

    if isinstance(data,pysh.shclasses.SHCoeffs):
        dataIn = XRcoeffsFromSH(data)
    else:
        dataIn = data.copy()

    dataIn = unstackSphDims(dataIn)

    pdTab, _ = multiDimXrToPD(dataIn.squeeze(drop=True).sortby(dataIn.attrs['harmonics']['dimList']), rowDims=dataIn.attrs['harmonics']['lDim'], fillna=fillna)

    return pdTab




def sphRealConvert(dataIn, method = 'sh', keyDims = None, incConj = True, rotPhase = None, addCSphase = False):
    """
    Convert real harmonics to complex form.

    Xarray for input and outputs.

    See also `PEMtk.sym.symHarm.setCoeffsSH`

    Parameters
    ----------
    dataIn : Xarray
        Xarray containing values to convert.
        Must contain keyDims (stacked or unstacked).
        TODO: generalise this, should get from :py:func:`epsproc.util.listFuncs.getRefDims()`.

    method : str, optional, default = 'sh'
        - 'sh': Set terms as per SHtools definitions, https://shtools.github.io/SHTOOLS/complex-spherical-harmonics.html
                Note this form has purely real or purely imaginary terms.
        - 'std': Set +/-m terms from |m| real harmonics as per wikipedia defn, https://en.wikipedia.org/wiki/Spherical_harmonics#Real_form
                Note: this form currently assumes only +M or -M for any given term, and also sets conjugates upon conversion if incConj=True.
                This may be phase-rotated relative to SHtools convention.
                NOTE: THIS METHOD CURRENTLY FAILS - phase issues somewhere, produces incorrect expansions in complex harmonics.

    keyDims : dict, optional, default = None
        Defines dim names for harmonic rank and order, and stacked dim name.
        If None, get standard names from YLMdimList() (=={'LM':['l','m']}).

    incConj : bool, optional, default = True
        Set opposite sign M term as conjugate in conversion.
        Generally should be True, but may cause issues in cases where both +/-M are defined in input.

    rotPhase : float, optional, default = None
        If not None, apply additional phase rotation by exp(-i*rotPhase), i.e. rotate phase about Z-axis.
        To match SHtools defn. set rotPhase = -np.pi/2 which will flip (x,y) axes.
        NOTE: CURRENTLY NOT CORRECT IN GENERAL

    addCSphase : bool, optional, default = False
        Apply additional Condon-Shortley (-1)^m term to swap phase/rotation direction?
        In general don't need this for consistency with Scipy derived harmonics, which already include CS phase.
        May need in some cases for consistency in (x,y) orientation of expanded distribution.


    Returns
    -------

    dataC : Xarray
        Results of conversion.
        All parameters set in dataC['harmonics']


    TODO

    - Generalise rotations, currently only supports theta (about Z-axis) to allow matching to SHtools conventions.
        See ep.sphCalc.TKQarrayRotX, which does work but needs some tidying up.
        E.g.    # Test existing frame rot... quite old code
                testBLMs = sphRealConvert(lmReal, method='std').rename({'LM':'BLM'})
                testBLMs.attrs['dataType'] = 'BLM'
                BLMrot, *_ = ep.TKQarrayRotX(testBLMs, ep.setPolGeoms())

    """

    dataCalc = dataIn.copy()
    # dataCalc.attrs = dataIn.attrs.copy()   # XR issue with attrs copy? Works for base dict, but not nested dicts.
    dataCalc.attrs = copy.deepcopy(dataIn.attrs)  # THIS WORKS ALSO FOR NESTED DICT CASE
                                                  # Definitely required in XR2022.3, 2022.6, but should be fixed in later versions, see https://github.com/pydata/xarray/issues/2835

    #*** Basic dim handling
    # TODO: should use checkDims here
    # LMStackFlag = False
    # stackDim = list(keyDims.keys())[0]
    # lDim = keyDims[stackDim][0]
    # mDim = keyDims[stackDim][1]
    #
    # if stackDim in dataCalc.dims:
    #     dataCalc =  dataCalc.unstack(stackDim)
    #     LMStackFlag = True

    # if keyDims is None:
    #     keyDims = YLMdimList(sType = 'sDict')   # Set to default dim names.

    # Functional version - returns above params to dataCalc.attrs['harmonics']
    # dataCalc = checkSphDims(dataCalc, keyDims)  # 22/09/22 - use separate check and unstack funs now
    dataCalc = unstackSphDims(dataCalc, keyDims)
    mDim = dataCalc.attrs['harmonics']['mDim']  # For brevity below.

    #*** METHODS
    if (method == 'std-v1') or (method == 'std'):
        # Set +/-m terms from |m| real harmonics as per
        # https://en.wikipedia.org/wiki/Spherical_harmonics#Real_form

        # Compute for all m, subselect +/- later.
        # plusMcase = ((-1)**np.abs(dataCalc.m))/np.sqrt(2)*(dataCalc + 1j*dataCalc)
        # plusMcase = ((-1)**np.abs(dataIn.m))/np.sqrt(2)*(dataIn + 1j*dataIn)  # For -ve m may need abs(m) here
        negMcase = ((-1)**np.abs(dataCalc[mDim]))/np.sqrt(2)*(dataCalc - 1j*dataCalc)
        # negMcase = (1/np.sqrt(2))*(dataCalc - 1j*dataCalc)
        # plusMcase = ((-1)**np.abs(dataCalc[mDim]))/np.sqrt(2)*(dataCalc + 1j*dataCalc)
        plusMcase = (1/np.sqrt(2))*(dataCalc + 1j*dataCalc)   # WITHOUT CS phase

        # Concat from various parts - OK.
        nonZeroM = xr.concat([plusMcase.where(plusMcase[mDim]>0, drop=True),  negMcase.where(negMcase[mDim]<0, drop=True)], dim=mDim)  # Keep all m in both cases

        # Set conjugate pairs?
        if incConj:
            dataC = xr.concat([dataCalc.sel({mDim:0}), nonZeroM, sphConj(nonZeroM)], dim=mDim)
        else:
            dataC = xr.concat([dataCalc.sel({mDim:0}), nonZeroM], dim=mDim)

    elif method == 'std-v2':
        # Set by coord assignment... have some norm issues with v1?
        # Test with subselection & different stacking...
        plusMterms = dataCalc.where(dataCalc[mDim]>0)
        negMterms = dataCalc.where(dataCalc[mDim]<0)

        plusMcase = ((-1)**np.abs(plusMterms[mDim]))/np.sqrt(2)*(plusMterms + 1j*plusMterms)
        # plusMcase = ((-1)**np.abs(dataIn.m))/np.sqrt(2)*(dataIn + 1j*dataIn)  # For -ve m may need abs(m) here
        negMcase = (1/np.sqrt(2))*(negMterms - 1j*negMterms)

        # Concat from various parts - OK.
        nonZeroM = xr.concat([plusMcase.where(plusMcase[mDim]>0, drop=True),  negMcase.where(negMcase[mDim]<0, drop=True)], dim=mDim)  # Keep all m in both cases

        # Set conjugate pairs?
        if incConj:
            dataC = xr.concat([dataCalc.sel({mDim:0}), nonZeroM, sphConj(nonZeroM)], dim=mDim)
        else:
            dataC = xr.concat([dataCalc.sel({mDim:0}), nonZeroM], dim=mDim)


    elif method == 'sh':
        # Set by coord assignment... Should match SHtools case.
        dataC = xr.zeros_like(dataCalc)
        # dataC = dataC.where(dataC[mDim]>-1, 1/np.sqrt(2)*(- 1j*dataCalc))   # -m case
        # dataC = dataC.where(dataC[mDim]>-1,((-1)**np.abs(dataCalc[mDim]))/np.sqrt(2)*(1j*dataCalc))   # -m case
        dataC = dataC.where(dataC[mDim]>-1,1/np.sqrt(2)*(1j*dataCalc))   # -m case, no additional (-1)^m phase term
        # dataC3 = dataC3.where(dataC3.m<1,((-1)**np.abs(dataIn.m))/np.sqrt(2)*(dataIn + 1j*dataIn))   # +m case
        dataC = dataC.where(dataC[mDim]<1,1/np.sqrt(2)*(dataCalc))   # +m case
        dataC = dataC.where(dataC[mDim]!=0,dataCalc)  # m=0 case

        if incConj:
            # Generate +/-m pairs from above
            # 27/09/22: added `stackOutput = False` to force dims to match input in updated sphConj() function.
            dataC = xr.concat([dataC, sphConj(dataC.where(dataC[mDim]!=0), stackOutput = False)], dim=mDim)  #, coords='minimal')
            # dataCp = dataCp.sortby(keyDims[stackDim]).stack(keyDims).dropna(dim=stackDim,how='all')


    else:
        print(f"*** sphRealConvert(method = {method}) not recognised.")

        return 0


    #*** Tidy up

    # Swap phase/rotation direction - need this for consistency with (x,y) orientation in current case.
    # May indicate inconsistently applied Condon-Shortley phase elsewhere?
    # UPDATE: yes, matches case for SHtools csphase = -1, which should be correct since this is included in Ylm definitions herein.
    # UPDATE2: Now fixed in real Ylm definitions.
    if addCSphase:
        dataC = (-1)**np.abs(dataC.m) * dataC   # Note abs(m) here, can't ** with -ve ints, but doesn't affect result in this case.

    # ADDITIONAL PHASE ROTATION
    if rotPhase is not None:
        # To match SHtools defn. set rotPhase = -np.pi/2 which will flip (x,y) axes.
        dataC = dataC * np.exp(-1j*rotPhase)

    # dataC = dataC.sortby(keyDims[stackDim]).stack(keyDims).dropna(dim=stackDim,how='all')  # Sort dims, stack and clean up.
    dataC = dataC.sortby(dataCalc.attrs['harmonics']['dimList']).stack(dataCalc.attrs['harmonics']['keyDims']).dropna(dim=dataCalc.attrs['harmonics']['stackDim'],how='all')  # Sort dims, stack and clean up.


    # Propagate attrs
    dataC.attrs.update(dataCalc.attrs)
    # dataC.attrs = dataCalc.attrs.copy()
    # # dataC = checkSphDims(dataC, keyDims)  # Will also unstack!
    dataC.attrs['harmonics'].update(YLMtype(method={'sphRealConvert':method},incConj=incConj))  #,keyDims=keyDims))

    # dataC.attrs['harmonics'] = {'dtype':'Complex harmonics',
    #                             'kind':'complex',
    #                             'method': method,
    #                             'incConj':incConj,
    #                             'keyDims':keyDims}

    return dataC


# DEV CODE
# def sphRealConvert(dataIn):
#     """
#     Convert real harmonics to complex form.
#
#     Xarray for input and outputs.
#
#     See also `PEMtk.sym.symHarm.setCoeffsSH`
#
#     """
#
#     # Set +/-m terms from |m| real harmonics as per
#     # https://en.wikipedia.org/wiki/Spherical_harmonics#Real_form
#
#     # Basic test - assume |m| > +/-m terms.
#     # May be incorrect phase for -m real harmonics?
#     # Note need abs(m), otherwise (-1)**m may throw errors
#     plusMcase = ((-1)**np.abs(dataIn.m))/np.sqrt(2)*(dataIn + 1j*dataIn)
#     # plusMcase = ((-1)**np.abs(dataIn.m))/np.sqrt(2)*(dataIn + 1j*dataIn)  # For -ve m may need abs(m) here
#     negMcase = (1/np.sqrt(2))*(dataIn - 1j*dataIn)
#     negMcase = negMcase.unstack('LM')
#     negMcase.coords['m'] = -negMcase.m
#     # negMcase = negMcase.stack({'LM':['l','m']})
#
#     # For case with ONLY + or -m in real harmonics (should be general...?) need to set by symmetry relations
#     # TODO: is this general? Should be since it's just a phase convention?
#     # THIS MATCHES WIKI defn., but in testing gives conjugate answer - might be accidental switch in numerics?
#     minM = dataIn.m.min()
#     if minM <0:
#         negMcase2 = (1/np.sqrt(2))*(dataIn - 1j*dataIn)
#         plusMcase2 = np.conj(((-1)**np.abs(dataIn.m)) *negMcase2)
#
#         negMcase2 = negMcase2.unstack('LM')
#         plusMcase2 = plusMcase2.unstack('LM')
#
#         plusMcase2.coords['m'] = -plusMcase2.m
#
#
#     if minM==0:
#         plusMcase2 = ((-1)**dataIn.m)/np.sqrt(2)*(dataIn + 1j*dataIn)
#
#         negMcase2 = np.conj(((-1)**np.abs(dataIn.m)) *plusMcase2)
#
#         plusMcase2 = plusMcase2.unstack('LM')
#         negMcase2 = negMcase2.unstack('LM')
#         negMcase2.coords['m'] = -negMcase2.m
#
#
#     # Concat - BASICALLY WORKING, just need to drop/fix m=0 terms.
#     # UPDATE - now fixed with .where()
#
#     # dataC = plusMcase + negMcase   # Keeps m=0 terms only!
#
#     # Concat from various parts - OK.
#     dataC = xr.concat([dataIn.unstack('LM').sel(m=0), plusMcase.where(plusMcase.m>0, drop=True).unstack('LM'),
#                        negMcase.where(negMcase.m<0, drop=True)], dim='m')
#     #Tidy up
#     dataC = dataC.stack({'LM':['l','m']}).dropna(dim='LM',how='all')
#     dataC.attrs['dtype'] = 'Complex harmonics'
#     dataC.attrs['kind'] = 'complex'
#     # dataC = dataC.stack({'LM':['l','m']}).dropna(dim='LM',how='all')
#
#
# #     return dataC, dataIn, plusMcase2, negMcase2
#
#     # Test 2nd form...
#     # Concat from various parts
#
#     # V1 - kills some terms due to removal of +/-m coords? May be able to fix with coords = something?
# #     dataC2 = xr.concat([dataIn.unstack('LM').sel(m=0), plusMcase2.where(plusMcase2.m>0, drop=True),
# #                        negMcase2.where(negMcase2.m<0, drop=True)], dim='m')  #, coords='minimal')
#
#     # V2 - correctly treats m=0 and keeps all other terms, not sure if other terms ARE CORRECT HOWEVER
#     dataC2 = xr.concat([dataIn.unstack('LM').sel(m=0),
#                     plusMcase2.where(plusMcase2.m != 0),
#                     negMcase2.where(negMcase2.m != 0)], dim='m')  #, coords='minimal')
#
#     #Tidy up
#     dataC2 = dataC2.stack({'LM':['l','m']}).dropna(dim='LM',how='all')
#     dataC2.attrs['dtype'] = 'Complex harmonics'
#     dataC2.attrs['kind'] = 'complex'
#
#     return dataC, dataC2
#
#     # V3 - swap on coords only
#     dataC3 = xr.zeros_like(dataIn)
#     dataC3.where(dataC3.m<0,(1/np.sqrt(2))*(dataIn - 1j*dataIn))
#
#
#
# #     # V4 - defn. from https://shtools.github.io/SHTOOLS/complex-spherical-harmonics.html
# #     # This defines coeff transforms, so should be correct.
# #     plusMcase = 1/np.sqrt(2)*(dataIn - 1j*dataIn)
# #     # plusMcase = ((-1)**np.abs(dataIn.m))/np.sqrt(2)*(dataIn + 1j*dataIn)  # For -ve m may need abs(m) here
# #     negMcase = (1/np.sqrt(2))*(dataIn - 1j*dataIn)
# #     negMcase = negMcase.unstack('LM')
# #     negMcase.coords['m'] = -negMcase.m
#



def sphConj(dataIn, stackOutput = True, cleanOutput = True, keepAttrs = True):
    """
    Compute conjugate harmonics from Xarray (includes +/-m sign switch).

    Yl,m = conj((-1)^m * Yl,-m)

    Basic dim handling with checkSphDims(dataIn) implemented, but should generalise.

    TODO: add switch for (-1)^m term - should drop for csphase=False case?

    """

    dataCalc = dataIn.copy()
    dataCalc.attrs = copy.deepcopy(dataIn.attrs)  # THIS WORKS ALSO FOR NESTED DICT CASE
                                                  # Definitely required in XR2022.3, 2022.6, but should be fixed in later versions, see https://github.com/pydata/xarray/issues/2835

    # LMStackFlag = False
    # if 'LM' in dataCalc.dims:
    #     dataCalc =  dataCalc.unstack('LM')
    #     LMStackFlag = True

    # Functional version - returns above params to dataCalc.attrs['harmonics']
    # Only run if not already set (e.g. already run in sphRealConvert(), including keyDims)
    # if (not 'harmonics' in dataCalc.attrs.keys()) or (not 'mDim' in dataCalc.attrs['harmonics'].keys()):
    #     dataCalc = checkSphDims(dataCalc)   # Just run with default keyDims - needs rationalisation.

    # Always run? This avoids issues with passed attrs in some cases (e.g. if set on a different dataIn style)
    # TODO: overwrite option?
    # dataCalc = checkSphDims(dataCalc)   # Just run with default keyDims - needs rationalisation.
    dataCalc = unstackSphDims(dataCalc)   # Not sure if required here? But might simplify restacking logic later.
                                          # UPDATE: now also runs checkSphDims if required.

    mDim = dataCalc.attrs['harmonics']['mDim']  # For brevity below.

    dataOut = np.conj((-1)**np.abs(dataCalc[mDim]) * dataCalc)    # Conj resultant
    # dataOut = (-1)**np.abs(dataCalc.m) * np.conj(dataCalc)  # Conj inputs only
    dataOut[mDim] = -1 * dataOut[mDim]

    if cleanOutput:
        dataOut = dataOut.dropna(dim=mDim,how='all')

    if stackOutput and dataCalc.attrs['harmonics']['LMStackFlag']:
        dataOut = dataOut.stack(dataCalc.attrs['harmonics']['keyDims'])

        # Tidy again?
        if cleanOutput:
            dataOut = dataOut.dropna(dim=dataCalc.attrs['harmonics']['stackDim'],how='all')

    if keepAttrs:
        dataOut.attrs = dataCalc.attrs

    return dataOut


def cleanLMcoords(daIn, refDims = None):
    """
    Basic routine to clean (L,M) coords to remove terms with m>l.

    Parameters
    ----------
    daIn : Xarray
        Data to clean.

    refDims : list, optional, default = None
        Default case checks uses checkSphDims(dataIn) to get (L,M) coords.
        Or or specify coord pair, e.g. ['l','m'].

    Returns
    -------
    daOut : Xarray with cleaned coords.
        Note in some cases this may not work, and stack/unstack cycles may reintroduce m>l terms.


    v1 29/03/23

    """

    da = daIn.copy()

    # Quick dim checks
    # For default case check for (l,m) and (L,M)
    # NOTE - in some updated ePSproc routines these are now set in da.attrs, should propagate this!
    if refDims is None:
        # Manual case
#         dimCheck = ep.util.misc.checkDims(da, refDims=['l','m'])

#         if not dimCheck['refDims']:
#             dimCheck = ep.util.misc.checkDims(da, refDims=['L','M'])

        # With function
        # da = ep.sphFuncs.sphConv.checkSphDims(da)
        da = checkSphDims(da)

        # Use da.attrs['harmonics']['dimList'] or 'lDim' and 'mDim'
        dimCheck = {'refDims': da.attrs['harmonics']['dimList']}

    else:
        dimCheck = ep.util.misc.checkDims(da, refDims=refDims)

    # Clean array - note this assumes [l,m] ordering.
    daOut = da.where(np.abs(da.coords[dimCheck['refDims'][1]])<=da.coords[dimCheck['refDims'][0]],drop=True)

    return daOut
