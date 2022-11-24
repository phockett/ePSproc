
import numpy as np
import xarray as xr # Currently used for type checks only.

# from epsproc.util import matEleSelector   # Circular/undefined import issue - call in function instead for now.
from epsproc import sphCalc
from epsproc.geomFunc import geomCalc
# from epsproc.geomFunc.geomCalc import (EPR, MFproj, betaTerm, remapllpL, w3jTable,)
from epsproc.geomFunc.geomUtils import genllpMatE, degenChecks

# Code as developed 16/17 March 2020.
# Needs some tidying, and should implement BLM Xarray attrs and format for output.
def afblmXprod(matEin, QNs = None, AKQS = None, EPRX = None, p=[0],
                BLMtable = None, BLMtableResort = None,
                lambdaTerm = None,
                # RX = None, eulerAngs = None, polLabel = None,
                polProd = None, AFterm = None,
                # basisDict = {},  May want to pass full dict here, or just pass as **basisDict from calling fn?
                thres = 1e-2, thresDims = 'Eke', selDims = {'Type':'L'}, sqThres = True, dropThres = True,  #, 'it':1},
                # sumDims = ['mu', 'mup', 'l','lp','m','mp'], sumDimsPol = ['P','R','Rp','p','S-Rp'], symSum = True,
                sumDims = ['mu', 'mup', 'l','lp','m','mp','S-Rp'], sumDimsPol = ['P','R','Rp','p'], symSum = True,  # Fixed summation ordering for AF*pol term...?
                degenDrop = True, SFflag = False, SFflagRenorm = False,
                BLMRenorm = 1,
                squeeze = False, phaseConvention = 'E',  #  , phaseCons = None
                basisReturn = "BLM", verbose = 0, **kwargs):
    r"""
    Implement :math:`\beta_{LM}^{AF}` calculation as product of tensors.

    .. math::
        \begin{eqnarray}
        \beta_{L,-M}^{\mu_{i},\mu_{f}} & =(-1)^{M} & \sum_{P,R',R}{[P]^{\frac{1}{2}}}{E_{P-R}(\hat{e};\mu_{0})}\sum_{l,m,\mu}\sum_{l',m',\mu'}(-1)^{(\mu'-\mu_{0})}{\Lambda_{R'}(\mu,P,R')B_{L,-M}(l,l',m,m')}I_{l,m,\mu}^{p_{i}\mu_{i},p_{f}\mu_{f}}(E)I_{l',m',\mu'}^{p_{i}\mu_{i},p_{f}\mu_{f}*}(E)\sum_{K,Q,S}\Delta_{L,M}(K,Q,S)A_{Q,S}^{K}(t)
        \end{eqnarray}


    Where each component is defined by fns. in :py:module:`epsproc.geomFunc.geomCalc` module.

    24/11/22 Added sqThres = True, dropThres = True for use with INITIAL matEleSelector() call only - in some cases may get issues with dim drops here otherwise.

    04/05/22 Added **kwargs, unused but allows for arb basis dict unpack and passing from other functions. May want to pipe back to Full basis return however.

    06/08/21 Added basic handling for degenerate states, including `degenDrop` option.
             Updated docs, but still rather messy!

    27/07/21 Removed eulerAngs & RX input options, since these are redundant (and lead to confusion here!).
             For cases where E-field and alignment distribution are rotated, set AKQS in rotated frame, see https://epsproc.readthedocs.io/en/latest/tests/ePSproc_frame_rotation_tests_Dec2019.html
             Also removed selDims={'it':1}, which can result in issues! In current code, adding 'it' to sumDims doesn't work (dim issues somewhere), so may be best to treat independently...?

    03/05/21 Tidying up a bit & improving/wrapping for fitting use (inc. basis function reuse).

    10/09/20 Verified (but messy) version, with updated defaults.

    01/09/20 Verified (but messy) version, including correct renormalisation routines.

    15/06/20 In progress!  Using mfblmXprod() as template, with just modified lambda term, and new alignment term, to change.


    For basic use see the docs: https://epsproc.readthedocs.io/en/dev/demos/ePSproc_class_demo_161020.html#Compute-LF/AF-\beta_{LM}-and-PADs


    Dev code:

        - afblmGeom_v1-ref_2020.py - Messy working v1 for reference, archived 04/05/21. Now working on tidier version.
        - geometric_method_dev_pt3_AFBLM_090620.ipynb
        - http://localhost:8888/lab/tree/dev/ePSproc/geometric_method_dev_Betas_090320.ipynb
        - D:\code\ePSproc\python_dev\ePSproc_MFBLM_Numba_dev_tests_120220.PY


    Parameters
    ----------
    matE : Xarray
        Xarray containing matrix elements, with QNs (l,m), as created by :py:func:`readMatEle`

    *** Optional calculation settings

    selDims : dict, default = {'Type':'L'}
        Selection parameters for calculations, may be be checked and appened herein.

    sumDims : list, default = ['mu', 'mup', 'l','lp','m','mp','S-Rp']
        Main summation dims, will be checked herein.

    sumDimsPol : list, default = ['P','R','Rp','p']
        Additional polarization summation dims.

    symSum : bool, default = True
        Sum over symmetries sets (=={Cont, Targ, Total}) if true.

    degenDrop : bool
        Flag to set dropping of degenerate components.

    thres : float, default = 1e-2
        Set threshold value, used to set input matrix elements and again for outputs.

    thresDims : str, default = 'Eke'
        Set threshold dimension (set to be contiguous).

    verbose : bool or int
        Print output?


    *** Optional renormalisation settings (mainly for testing only)

    SFflag : bool, default = False
        Multiply input matrix elements by complex scale-factor if true.

    SFflagRenorm : bool, default = False
        Renorm output BLMs by complex scale-factor if true.

    BLMRenorm : int, default = 1
        Set different BLM renorm conventions.
        If 1 renorm by B00.
        See  code for further details.

    sqThres  : bool, default = True
        Used for initial matrix element threholding call only.
        Default = True as previous code, but can cause issues in custom cases (usually accidentally dropping singleton mu coord).
        Added 24/11/22

    dropThres  : bool, default = True
        Used for initial matrix element threholding call only.
        Default = True as previous code, but can cause issues in custom cases (usually accidentally dropping singleton mu coord).
        Added 24/11/22

    squeeze : bool, default = False
        Squeeze output array after thresholding?
        Note: this may cause dim issues if True.


    *** Optional input data/basis functions (mainly for fitting routine use)

    QNs : np.array, optional, default = None
        List of QNs as generated by :py:func:`genllpMatE`.
        Will be generated if not passed.

    AKQS : Xarray, optional, default = None
        Alignment parameters, as set by :py:func:`setADMs`.
        Defaults to isotropic case if not passed.

    EPRX : Xarray, optional, default = None
        E-field parameters, as generated by :py:func:`EPR`.
        Defaults to normalised/unity field, pol = p (below).

    p : list or array, optional, default = [0]
        Specify polarization terms p.
        Possibly currently only valid for p=0, TBC
        See https://epsproc.readthedocs.io/en/latest/methods/geometric_method_dev_260220_090420_tidy.html#E_{P,R}-tensor

    BLMtable, BLMtableResort : Xarrays, optional, default = None
        Beta calculation parameters, as defined by :py:func:`geomCalc.betaTerm`.
        BLMtableResort includes phase settings & param renaming as herein.

    lambdaTerm : Xarray, optional, default = None
        Lambda term parameters, as defined by :py:func:`geomCalc.MFproj`

    AFterm : Xarray, optional, default = None
        Alignment term as defined by :py:func:`geomCalc.deltaLMKQS`

    polProd : Xarray, optional, default = None
        Polarization tensor as defined by EPRXresort * lambdaTermResort * AFterm

    phaseConvention : optional, str, default = 'E'
        Set phase conventions with :py:func:`epsproc.geomCalc.setPhaseConventions`.
        To use preset phase conventions, pass existing dictionary.

    basisReturn : optional, str, default = "BLM"
        - 'BLM' return Xarray of results only.
        - 'Full' return Xarray of results + basis set dictionary as set during the run.
        - 'Product', as full, but minimal basis set with products only.
        - 'Results' or 'Legacy' direct return of various calc. results Xarrays.

    **kwargs, unused but allows for arb basis dict unpack and passing from other functions.

    Returns
    -------
    Xarray
        Set of AFBLM calculation results

    dict
        Dictionary of basis functions, only if basisReturn != 'BLM' (see basisReturn paramter notes).


    Notes
    -----

    Cross-section outputs now set as:

    - XSraw = direct AF calculation output.
    - XSrescaled = XSraw * sqrt(4pi)
    - XSiso = direct sum over matrix elements

    Where XSrescaled == XSiso == ePS GetCro output for isotropic distribution.
    Optionally set SFflag = True to multiply by (complex) scale-factor.


    """
    from epsproc.util import matEleSelector

    calcSettings = locals()  # Grab passed args (calc. settings) for reference later.

    # Set phase conventions - either from function call or via passed dict.
    # if type(phaseConvention) is str:
    #     phaseCons = geomCalc.setPhaseConventions(phaseConvention = phaseConvention)
    # else:
    #     phaseCons = phaseConvention

    # For transparency/consistency with subfunctions, str/dict now set in setPhaseConventions()
    phaseCons = geomCalc.setPhaseConventions(phaseConvention = phaseConvention)

    # Fudge - set this for now to enforce additonal unstack and phase corrections later.
    # 03/05/21 - now in passed args for basis set.
    # BLMtableResort = None

    #*** Threshold and selection
    # Make explicit copy of data to avoid any overwrite issues
    matE = matEin.copy()
    matE.attrs = matEin.attrs  # May not be necessary with updated Xarray versions

    # Use SF (scale factor)
    # Write to data.values to make sure attribs are maintained. (Not the case for da = da*da.SF)
    if SFflag:
        matE.values = matE * matE.SF

    # Degenerate state handling
    degenDict = degenChecks(matE, selDims, sumDims, degenDrop, verbose)
    selDims = degenDict['selDims']  # Update selDims
    calcSettings['degenDict'] = degenDict  # Add to calcSettings for output with main Xarray attribs.

    # Threshold
    matEthres = matEleSelector(matE, thres = thres, inds = selDims, dims = thresDims, sq = sqThres, drop = dropThres)

    # Sum **AFTER** threshold and selection, to allow for subselection on symmetries via matEleSelector
    symDegen = 1
    if 'Sym' in matEthres.dims:
        symDegen = matEthres.Sym.size  # Set degeneracy - use thresholded or raw matrix elements here?

        if symSum:
            matEthres = matEthres.sum('Sym')  # Sum over ['Cont','Targ','Total'] stacked dims.


    #*** Polarization terms
    if (EPRX is None) and (polProd is None):  # Skip if product term already passed
        # *** EPR
        # EPRX = geomCalc.EPR(form = 'xarray', p = p, phaseConvention = phaseConvention).sel({'R-p':0})  # Set for R-p = 0 for p=0 case (redundant coord) - need to fix in e-field mult term!
        # EPRXresort = EPRX.unstack().squeeze().drop('l').drop('lp')  # This removes photon (l,lp) dims fully. Be careful with squeeze() - sends singleton dims to non-dimensional labels.
#         EPRXresort = EPRX.unstack().drop('l').drop('lp')  # This removes photon (l,lp) dims fully, but keeps (p,R) as singleton dims.
#         EPRXresort = EPRX.unstack().squeeze(['l','lp']).drop(['l','lp'])  # Safe squeeze & drop of selected singleton dims only.

#         EPRX = geomCalc.EPR(form = 'xarray', p = p).unstack().sum(['p','R-p'])  # Set for general sum over (p,R-p) terms - STILL need to fix in e-field mult term!
#         EPRX = geomCalc.EPR(form = 'xarray', p = p).unstack().sum('R-p')  # Set for general sum over (p,R-p) terms - STILL need to fix in e-field mult term!
# TODO: check and fix if necessary for p!=0 case
        EPRX = geomCalc.EPR(form = 'xarray', p = p, phaseConvention = phaseConvention).unstack().sel({'R-p':0}).drop('R-p')  # Working case as of v1.3.0-dev, but valid for p=0 only?
        EPRXresort = EPRX.squeeze(['l','lp']).drop(['l','lp'])  # Safe squeeze & drop of selected singleton dims only.

        if phaseCons['mfblmCons']['negRcoordSwap']:
            EPRXresort['R'] *= -1

    if (lambdaTerm is None) and (polProd is None):  # Skip if product term already passed
        # Set polGeoms if Euler angles are passed.
        # if eulerAngs is not None:

        # 27/07/21 - removed extraneous (and possibly erroneous) frame rotation "option"
        # Set explictly here - only want (0,0,0) term in any case!
        # eulerAngs = np.array([0,0,0], ndmin=2)
        # RX = ep.setPolGeoms(eulerAngs = eulerAngs)   # This throws error in geomCalc.MFproj???? Something to do with form of terms passed to wD, line 970 vs. 976 in geomCalc.py
        # Set polGeoms if Euler angles are passed.
        # if eulerAngs is not None:
        #     RX = setPolGeoms(eulerAngs = eulerAngs)
        #
        # if RX is None:
        #     # Alternatively - just set default values then sub-select.
        #     RX = sphCalc.setPolGeoms()
        #
        # Subselect on pol geoms if label is passed.
        # May want to add try/fail here as this might be a bit fragile.
        # if polLabel is not None:
        #     RX = RX.sel({'Label':polLabel})

        RX = sphCalc.setPolGeoms(eulerAngs = [0,0,0])  # Use setPolGeoms, but ONLY VALID FOR (0,0,0) case BY DEFINITION (no frame rotation term in AF formulation, although can ACCIDENTALLY APPLY with MFproj() function below).

        # *** Lambda term
        lambdaTerm, lambdaTable, lambdaD, _ = geomCalc.MFproj(RX = RX, form = 'xarray', phaseConvention = phaseConvention)
        # lambdaTermResort = lambdaTerm.squeeze().drop('l').drop('lp')   # This removes photon (l,lp) dims fully.
        # lambdaTermResort = lambdaTerm.squeeze(['l','lp']).drop(['l','lp'])  # Safe squeeze & drop of selected singleton dims only.
        # lambdaTermResort = lambdaTerm.squeeze(['l','lp']).drop(['l','lp']).sel({'Labels':'z'}).sum('R')  # Safe squeeze & drop of selected singleton dims only, select (0,0,0) term only for pol. geometry.
        lambdaTermResort = lambdaTerm.squeeze(['l','lp']).drop(['l','lp']).sum('R') # Without explicit geom selection.
        # NOTE dropping of redundant R coord here - otherwise get accidental R=Rp correlations later!

    # *** Blm term with specified QNs
    if (BLMtable is None) and (BLMtableResort is None):  # Skip this is BLMtableResort is passed

        # Set terms if not passed to function
        if QNs is None:
            QNs = genllpMatE(matEthres, phaseConvention = phaseConvention)

        QNsBLMtable = QNs.copy()

        # Switch signs (m,M) before 3j calcs.
        if phaseCons['mfblmCons']['BLMmPhase']:
            QNsBLMtable[:,3] *= -1
            QNsBLMtable[:,5] *= -1

        # QNsBLMtable[:,3] *= -1
        BLMtable = geomCalc.betaTerm(QNs = QNsBLMtable, form = 'xdaLM', phaseConvention = phaseConvention)

#         if BLMmPhase:
#             BLMtable['m'] *= -1

    if BLMtableResort is None:
        # Apply additional phase convention
        BLMtableResort = BLMtable.copy().unstack()

        if phaseCons['mfblmCons']['negMcoordSwap']:
            BLMtableResort['M'] *= -1

        if phaseCons['mfblmCons']['Mphase']:
            BLMtableResort *= np.power(-1, np.abs(BLMtableResort.M))  # Associated phase term

        if phaseCons['mfblmCons']['negmCoordSwap']:
            BLMtableResort['m'] *= -1

        if phaseCons['mfblmCons']['mPhase']:
            BLMtableResort *= np.power(-1, np.abs(BLMtableResort.m))  # Associated phase term

        # RENAME, M > (S-R') for AF case - this correctly allows for all MF projections!!!
        # Some MF phase cons as applied above may also be incorrect?
        BLMtableResort = BLMtableResort.rename({'M':'S-Rp'})


    #*** Alignment term
    if (AFterm is None) and (polProd is None):  # Skip if product term already passed

        if AKQS is None:
            AKQS = sphCalc.setADMs()     # If not passed, set to defaults - A(0,0,0)=1 term only, i.e. isotropic distribution.

        AFterm, DeltaKQS = geomCalc.deltaLMKQS(EPRXresort, AKQS, phaseConvention = phaseConvention)


    #*** Products

    # polProd, takes account of polarization + alignment (geometric) terms inc. sum over `sumDimsPol`
    if polProd is None:
        polProd = (EPRXresort * lambdaTermResort * AFterm)

        # Set additional phase term, (-1)^(mup-p) **** THIS MIGHT BE SPURIOUS FOR GENERAL EPR TENSOR CASE??? Not sure... but definitely won't work if p summed over above!
        if phaseCons['mfblmCons']['mupPhase']:
            mupPhaseTerm = np.power(-1, np.abs(polProd.mup - polProd.p))
            polProd *= mupPhaseTerm

        # Additional [P]^1/2 degen term, NOT included in EPR defn.
        # Added 09/04/20
        polProd *= np.sqrt(2*polProd.P+1)

        polProd = polProd.sum(sumDimsPol)

        polProd = matEleSelector(polProd, thres = thres)  # Select over dims for reduction.


    # Matrix element pair-wise multiplication by (l,m,mu) dims
    matEconj = matEthres.copy().conj()
    # matEconj = matEconj.unstack().rename({'l':'lp','m':'mp','mu':'mup'})  # Full unstack
    # matEmult = matEthres.unstack() * matEconj
    matEconj = matEconj.unstack('LM').rename({'l':'lp','m':'mp','mu':'mup'})  # Unstack LM only.
    matEmult = matEthres.unstack('LM') * matEconj
    matEmult.attrs['dataType'] = 'multTest'

    # Threshold product and drop dims.
    # matEmult = ep.util.matEleSelector(matEmult, thres = thres, dims = thresDims)
    matEmult = matEleSelector(matEmult, thres = thres, dims = thresDims)

    # Apply additional phase conventions?
    if phaseCons['afblmCons']['llpPhase']:
        matEmult *= np.power(-1, np.abs(matEmult.l - matEmult.lp))

    # Product terms with similar dims
    BLMprod = matEmult * BLMtableResort  # Unstacked case with phase correction - THIS DROPS SYM TERMS? Takes intersection of das - http://xarray.pydata.org/en/stable/computation.html#automatic-alignment


    # Test big mult...
    # mTerm = polProd.sel({'R':0,'Labels':'z'}) * BLMprod.sum(['Total'])    # With selection of z geom.  # BLMprod.sum(['Cont', 'Targ', 'Total'])
    # mTerm = polProd.sel({'R':0}) * BLMprod    # BLMprod.sum(['Cont', 'Targ', 'Total'])
    mTerm = polProd * BLMprod
    # Multiplication works OK, and is fast... but might be an ugly result... INDEED - result large and slow to manipulate, lots of dims and NaNs. Better to sub-select terms first!



    # No subselection, mTerm.size = 6804000
    # For polProd.sel({'R':0}), mTerm.size = 1360800
    # For polProd.sel({'R':0,'Labels':'z'}), mTerm.size = 453600
    # Adding also BLMprod.sum(['Total']), mTerm.size = 226800
    # Adding also BLMprod.sum(['Cont', 'Targ', 'Total']), mTerm.size = 113400  So, for sym specific calcs, may be better to do split-apply type methods


    # mTerm.attrs['file'] = 'MulTest'  # Temporarily adding this, not sure why this is an issue here however (not an issue for other cases...)
    mTerm.attrs = matEin.attrs  # Propagate attrs from input matrix elements.
    # mTerm.attrs['phaseConvention'] = {phaseConvention:phaseCons}  # Log phase conventions used.
    mTerm.attrs['phaseCons'] = geomCalc.setPhaseConventions(phaseConvention = phaseConvention)  # Log phase conventions used.


    # Sum and threshold
#     sumDims = ['P', 'mu', 'mup', 'Rp', ]  # Define dims to sum over
    xDim = {'LM':['L','M']}
    mTermSum = mTerm.sum(sumDims)

    if squeeze is True:
        mTermSum = mTermSum.squeeze()  # Leave this as optional, since it can cause issues for M=0 only case

    mTermSumThres = matEleSelector(mTermSum.stack(xDim), thres=thres, dims = thresDims)
#     mTermSumThres = mTermSum

    #*** Normalise

    # Additional factors & renorm - calc. XS as per lfblmGeom.py testing, verified vs. ePS outputs for B2 case, June 2020
    #  XSmatE = (matE * matE.conj()).sel(selDims).sum(['LM','mu'])  # (['LM','mu','it'])  # Cross section as sum over mat E elements squared (diagonal terms only)
    XSmatE = (matEthres * matEthres.conj()).sum(['LM','mu']) # .expand_dims({'t':[0]})  # Use selected & thresholded matE.
                                # NOTE - this may fail in current form if dims are missing.
                                # Quick hack for testing, add expand_dims({'t':[0]}) need to ensure matching dims for division!
    normBeta = 3/5 * (1/XSmatE)  # Normalise by sum over matrix elements squared.

    # Additional scaling if required for degeneracy and/or SF
    if degenDict['degenFlag']:
        mTermSumThres.values = mTermSumThres * degenDict['degenN']

    if SFflagRenorm:
        mTermSumThres.values = mTermSumThres/mTermSumThres.SF

    mTermSumThres['XSraw'] = mTermSumThres.sel({'L':0,'M':0}).drop('LM').copy()  # This basically works, and keeps all non-summed dims... but may give issues later...? Make sure to .copy(), otherwise it's just a pointer.

    # Rescale by sqrt(4pi)*SF, this matches GetCro XS outputs in testing.
    # mTermSumThres['XSrescaled'] = mTermSumThres['XSraw']*mTermSumThres['SF']*np.sqrt(4*np.pi)
    mTermSumThres['XSrescaled'] = mTermSumThres['XSraw']*np.sqrt(4*np.pi)

    # In some cases may also need to account for degen...? Seemed to in N2 AF testing 10/09/20, but may have been spurious result.
    # Could also be Sph <> Lg conversion issue?
    # if symSum:
    #     # Rescale by sqrt(4pi)*SF, this matches GetCro XS outputs in testing.
    #     mTermSumThres['XSrescaled'] = mTermSumThres['XSraw']*mTermSumThres['SF']*np.sqrt(4*np.pi)
    #
    # else:
    #     # mTermSumThres['XSrescaled'] /= symDegen  # Correct sym unsummed case (multiple summation issue?)
    #     # Actually, looks like issue is scaling for SF - for single sym case DON'T NEED IT to match GetCro outputs.
    #     # Is this then correct?
    #     mTermSumThres['XSrescaled'] = mTermSumThres['XSraw']*np.sqrt(4*np.pi)


    mTermSumThres['XSiso'] = XSmatE/3  # ePolyScat defn. for LF cross-section. (i.e. isotropic distribution)
    # mTermSumThres['XS2'] = symDegen * XSmatE/3  # Quick hack for testing, with symDegen

    # Renorm betas by B00?
    if BLMRenorm:
        # mTermSumThres /= mTermSumThres.sel({'L':0,'M':0}).drop('LM')
        if BLMRenorm == 0:
            # Keep values, scale by normBeta & sqrt(4pi), to match ePS values.
            mTermSumThres *= normBeta * np.sqrt(4*np.pi)

        if BLMRenorm == 1:
            # Renorm by full t-dependent XS only
            mTermSumThres /= mTermSumThres['XSraw']

        elif BLMRenorm == 2:
            # Renorm by isotropic XS only
            mTermSumThres /= mTermSumThres['XSiso']

        elif BLMRenorm == 3:
            # Renorm by isotropic XS, then t-dependent (calculated) XS, then additional factors.
            # mTermSumThres /= mTermSumThres['XSiso']  # Includes 1/3 norm factor
            mTermSumThres /= XSmatE
            mTermSumThres['XSrenorm'] = mTermSumThres.sel({'L':0,'M':0}).drop('LM').copy() # Enforce dims here, otherwise get stray L,M coords.
            mTermSumThres /= mTermSumThres['XSrenorm']
            # mTermSumThres *= symDegen/(2*mTermSumThres.L + 1)  # Renorm to match ePS GetCro defns. Not totally sure if symDegen is correct - TBC.
            # mTermSumThres *= symDegen/5  # Check if 2L+1 factor is correct...? This seems better for N2 AF test case, otherwise L>2 terms very small - maybe M-state degen only by matrix elements?
            # mTermSumThres /= (2*mTermSumThres.L + 1)
            # mTermSumThres = symDegen/5 * mTermSumThres.where(mTermSumThres.L > 0)
            # mTermSumThres = mTermSumThres.where(mTermSumThres.L == 0, symDegen/5 * mTermSumThres)
            # mTermSumThres = mTermSumThres.where(mTermSumThres.L == 0, symDegen/(2*mTermSumThres.L + 1) * mTermSumThres)
            mTermSumThres *= symDegen/(2*mTermSumThres.L + 1)

        elif BLMRenorm == 4:
            # Alt scheme... similar to #3, but testing different renorm factors
            # mTermSumThres /= mTermSumThres['XSiso'] # Includes 1/3 norm factor
            mTermSumThres /= XSmatE
            mTermSumThres['XSrenorm'] = mTermSumThres.sel({'L':0,'M':0}).drop('LM').copy() # Enforce dims here, otherwise get stray L,M coords.
            mTermSumThres /= mTermSumThres['XSrenorm']

            mTermSumThres *= symDegen
            mTermSumThres /= (2*mTermSumThres.L + 1)

        else:
            mTermSumThres *= normBeta  # Scale by normBeta only.

    # Propagate attrs
    mTermSum.attrs = mTerm.attrs
    mTermSum.attrs['dataType'] = 'multTest'
    mTermSum.attrs['BLMRenorm'] = BLMRenorm

    mTermSumThres.attrs = mTerm.attrs
    mTermSumThres.attrs['dataType'] = 'multTest'
    mTermSum.attrs['BLMRenorm'] = BLMRenorm

# TODO: Set XS as per old mfpad()
#     BLMXout['XS'] = (('Eke','Euler'), BLMXout[0].data)  # Set XS = B00
#     BLMXout = BLMXout/BLMXout.XS  # Normalise

    #**** Tidy up and reformat to standard BLM array (see ep.util.BLMdimList() )
    # TODO: finish this, and set this as standard output
    BetasNormX = mTermSumThres.unstack().rename({'L':'l','M':'m'}).stack({'BLM':['l','m']})

    # Set/propagate global properties
    BetasNormX.attrs = matE.attrs
    # BetasNormX.attrs['thres'] = thres

    # TODO: update this for **vargs
    # BLMXout.attrs['sumDims'] = sumDims # May want to explicitly propagate symmetries here...?
    # BLMXout.attrs['selDims'] = [(k,v) for k,v in selDims.items()]  # Can't use Xarray to_netcdf with dict set here, at least for netCDF3 defaults.

    # 28/07/21 added locals() to pipe all args > attrs. Ignore dict issue for now!
    # BetasNormX.attrs.update(calcSettings)  # This works, but can be a bit of a mess due to passed basis sets
    [BetasNormX.attrs.update({k:calcSettings[k]}) for k in calcSettings.keys() if not (isinstance(calcSettings[k], xr.DataArray))]  # Slightly ugly, but set only items which are not Xarrays.


    BetasNormX.attrs['dataType'] = 'BLM'

    # Set return args based on basisReturn parameter
    # Full results set, including all versions
    if verbose:
        print(f"Return type {basisReturn}.")

    if basisReturn in ["Results", "Legacy"]:
        # print("Legacy")
        return mTermSumThres, mTermSum, mTerm, BetasNormX

    # Return basis arrays/tensors
    elif basisReturn == "Full":
        # print("Full")
        basis = {'QNs':QNs, 'EPRX':EPRXresort, 'lambdaTerm':lambdaTermResort,
                'BLMtable':BLMtable, 'BLMtableResort':BLMtableResort,
                'AFterm':AFterm, 'AKQS':AKQS, 'polProd':polProd,
                'phaseConvention':phaseCons, 'BLMRenorm':BLMRenorm}   #, 'phaseCons':phaseCons}

        return  BetasNormX, basis

    # Return product basis fns. for use in fitting routines
    elif basisReturn == "ProductBasis":
        basis = {'BLMtableResort':BLMtableResort, 'polProd':polProd, 'phaseConvention':phaseCons, 'BLMRenorm':BLMRenorm}

        return  BetasNormX, basis

    # Minimal return
    elif basisReturn == "BLM":
        # print("BLM")
        return BetasNormX

    else:
        print(f"Return type {basisReturn} not recognised, defaulting to BLM only.")
        return BetasNormX


def AFwfExp(matE, AKQS=None, thres = 1e-2, thresDims = 'Eke', selDims = {'Type':'L', 'it':1}):
    r"""
    Implement (approximate) LF/AF wavefunction expansion,

    .. math::
        \begin{equation}
        ^{AF}T_{\mu_{0}}^{p_{i}\mu_{i},p_{f}\mu_{f}}(\hat{k}_{L})=8\pi^{2}\sum_{K,Q,S}\sum_{l,m,\mu,\Lambda}A_{Q,S}^{K}I_{l,m,\mu}^{p_{i}\mu_{i},p_{f}\mu_{f}}(E)(-1)^{m-\Lambda}(-1)^{\mu-\mu_{0}}(-1)^{Q-S}\left(\begin{array}{ccc}
        l & 1 & K\\
        -m & -\mu & -Q
        \end{array}\right)\left(\begin{array}{ccc}
        l & 1 & K\\
        -\Lambda & -\mu_{0} & -S
        \end{array}\right)Y_{l\Lambda}^{*}(\hat{k}_{M})
        \end{equation}

    Where each component is defined by fns. in :py:module:`epsproc.geomFunc.geomCalc` module.

    01/02/21 version 1, this is pretty rough and ready, and possibly not correct.

    See http://localhost:8888/lab/tree/SLAC/angular_streaking/AF_wavefns_method_dev_050121.ipynb for dev notes.

    """

    from epsproc.util import matEleSelector

    matEthres = matEleSelector(matE, thres = thres, inds = selDims, dims = thresDims, sq = True, drop = True)

    #*** Alignment term
    if AKQS is None:
        AKQS = sphCalc.setADMs()    # If not passed, set to defaults - A(0,0,0)=1 term only, i.e. isotropic distribution.

    # New alignment function
    AFterm, DeltaKQS = geomCalc.deltalmKQSwf(matEthres, AKQS) #, phaseConvention = phaseConvention)

    AFexp = matEthres.unstack('LM') * AFterm
    AFexp.attrs['dataType'] = 'AFwfExp'  # Set this for lmPlot

    return AFexp, AFterm, DeltaKQS # Return other values for debug
