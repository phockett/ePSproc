
import numpy as np

# from epsproc.util import matEleSelector   # Circular/undefined import issue - call in function instead for now.
from epsproc import sphCalc
from epsproc.geomFunc import geomCalc
# from epsproc.geomFunc.geomCalc import (EPR, MFproj, betaTerm, remapllpL, w3jTable,)
from epsproc.geomFunc.geomUtils import genllpMatE

# Code as developed 16/17 March 2020.
# Needs some tidying, and should implement BLM Xarray attrs and format for output.
def afblmXprod(matEin, QNs = None, AKQS = None, EPRX = None, p=[0], BLMtable = None,
                lambdaTerm = None, RX = None, eulerAngs = None,
                thres = 1e-2, thresDims = 'Eke', selDims = {'it':1, 'Type':'L'},
                # sumDims = ['mu', 'mup', 'l','lp','m','mp'], sumDimsPol = ['P','R','Rp','p','S-Rp'], symSum = True,
                sumDims = ['mu', 'mup', 'l','lp','m','mp','S-Rp'], sumDimsPol = ['P','R','Rp','p'], symSum = True,  # Fixed summation ordering for AF*pol term...?
                SFflag = False, SFflagRenorm = False, BLMRenorm = 1,
                squeeze = False, phaseConvention = 'S'):
    r"""
    Implement :math:`\beta_{LM}^{AF}` calculation as product of tensors.

    .. math::
        \begin{eqnarray}
        \beta_{L,-M}^{\mu_{i},\mu_{f}} & =(-1)^{M} & \sum_{P,R',R}{[P]^{\frac{1}{2}}}{E_{P-R}(\hat{e};\mu_{0})}\sum_{l,m,\mu}\sum_{l',m',\mu'}(-1)^{(\mu'-\mu_{0})}{\Lambda_{R'}(\mu,P,R')B_{L,-M}(l,l',m,m')}I_{l,m,\mu}^{p_{i}\mu_{i},p_{f}\mu_{f}}(E)I_{l',m',\mu'}^{p_{i}\mu_{i},p_{f}\mu_{f}*}(E)\sum_{K,Q,S}\Delta_{L,M}(K,Q,S)A_{Q,S}^{K}(t)
        \end{eqnarray}


    Where each component is defined by fns. in :py:module:`epsproc.geomFunc.geomCalc` module.

    10/09/20 Verified (but messy) version, with updated defaults.

    01/09/20 Verified (but messy) version, including correct renormalisation routines.

    15/06/20 In progress!  Using mfblmXprod() as template, with just modified lambda term, and new alignment term, to change.

    Dev code:
        geometric_method_dev_pt3_AFBLM_090620.ipynb
        http://localhost:8888/lab/tree/dev/ePSproc/geometric_method_dev_Betas_090320.ipynb
        D:\code\ePSproc\python_dev\ePSproc_MFBLM_Numba_dev_tests_120220.PY
    TOTAL MESS AT THE MOMENT>>?>>>>?DFdas<>r ty



    Parameters
    ----------

    phaseConvention : optional, str, default = 'S'
        Set phase conventions with :py:func:`epsproc.geomCalc.setPhaseConventions`.
        To use preset phase conventions, pass existing dictionary.



    Notes
    -----

    Cross-section outputs now set as:

    - XSraw = direct AF calculation output.
    - XSrescaled = XSraw * SF * sqrt(4pi)
    - XSiso = direct sum over matrix elements

    Where XSrescaled == XSiso == ePS GetCro output for isotropic distribution.

    """
    from epsproc.util import matEleSelector

    # Set phase conventions - either from function call or via passed dict.
    # if type(phaseConvention) is str:
    #     phaseCons = geomCalc.setPhaseConventions(phaseConvention = phaseConvention)
    # else:
    #     phaseCons = phaseConvention

    # For transparency/consistency with subfunctions, str/dict now set in setPhaseConventions()
    phaseCons = geomCalc.setPhaseConventions(phaseConvention = phaseConvention)

    # Fudge - set this for now to enforce additonal unstack and phase corrections later.
    BLMtableResort = None

    #*** Threshold and selection
    # Make explicit copy of data to avoid any overwrite issues
    matE = matEin.copy()
    matE.attrs = matEin.attrs  # May not be necessary with updated Xarray versions

    # Use SF (scale factor)
    # Write to data.values to make sure attribs are maintained. (Not the case for da = da*da.SF)
    if SFflag:
        matE.values = matE * matE.SF

    matEthres = matEleSelector(matE, thres = thres, inds = selDims, dims = thresDims, sq = True, drop = True)

    # Sum **AFTER** threshold and selection, to allow for subselection on symmetries via matEleSelector
    symDegen = 1
    if 'Sym' in matEthres.dims:
        symDegen = matEthres.Sym.size  # Set degeneracy - use thresholded or raw matrix elements here?

        if symSum:
            matEthres = matEthres.sum('Sym')  # Sum over ['Cont','Targ','Total'] stacked dims.


    # Set terms if not passed to function
    if QNs is None:
        QNs = genllpMatE(matEthres, phaseConvention = phaseConvention)

    #*** Polarization terms
    if EPRX is None:
        # *** EPR
        # EPRX = geomCalc.EPR(form = 'xarray', p = p, phaseConvention = phaseConvention).sel({'R-p':0})  # Set for R-p = 0 for p=0 case (redundant coord) - need to fix in e-field mult term!
        # EPRXresort = EPRX.unstack().squeeze().drop('l').drop('lp')  # This removes photon (l,lp) dims fully. Be careful with squeeze() - sends singleton dims to non-dimensional labels.
#         EPRXresort = EPRX.unstack().drop('l').drop('lp')  # This removes photon (l,lp) dims fully, but keeps (p,R) as singleton dims.
#         EPRXresort = EPRX.unstack().squeeze(['l','lp']).drop(['l','lp'])  # Safe squeeze & drop of selected singleton dims only.

#         EPRX = geomCalc.EPR(form = 'xarray', p = p).unstack().sum(['p','R-p'])  # Set for general sum over (p,R-p) terms - STILL need to fix in e-field mult term!
#         EPRX = geomCalc.EPR(form = 'xarray', p = p).unstack().sum('R-p')  # Set for general sum over (p,R-p) terms - STILL need to fix in e-field mult term!
        EPRX = geomCalc.EPR(form = 'xarray', p = p).unstack().sel({'R-p':0}).drop('R-p')
        EPRXresort = EPRX.squeeze(['l','lp']).drop(['l','lp'])  # Safe squeeze & drop of selected singleton dims only.

        if phaseCons['mfblmCons']['negRcoordSwap']:
            EPRXresort['R'] *= -1

    if lambdaTerm is None:
        # Set polGeoms if Euler angles are passed.
        # if eulerAngs is not None:

        # Set explictly here - only want (0,0,0) term in any case!
        # eulerAngs = np.array([0,0,0], ndmin=2)
        # RX = ep.setPolGeoms(eulerAngs = eulerAngs)   # This throws error in geomCalc.MFproj???? Something to do with form of terms passed to wD, line 970 vs. 976 in geomCalc.py
        # Alternatively - just set default values then sub-select.
        RX = sphCalc.setPolGeoms()

        # *** Lambda term
        lambdaTerm, lambdaTable, lambdaD, _ = geomCalc.MFproj(RX = RX, form = 'xarray', phaseConvention = phaseConvention)
        # lambdaTermResort = lambdaTerm.squeeze().drop('l').drop('lp')   # This removes photon (l,lp) dims fully.
        # lambdaTermResort = lambdaTerm.squeeze(['l','lp']).drop(['l','lp'])  # Safe squeeze & drop of selected singleton dims only.
        lambdaTermResort = lambdaTerm.squeeze(['l','lp']).drop(['l','lp']).sel({'Labels':'z'}).sum('R')  # Safe squeeze & drop of selected singleton dims only, select (0,0,0) term only for pol. geometry.
        # NOTE dropping of redundant R coord here - otherwise get accidental R=Rp correlations later!

    # *** Blm term with specified QNs
    if BLMtable is None:

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
    if AKQS is None:
        AKQS = sphCalc.setADMs()     # If not passed, set to defaults - A(0,0,0)=1 term only, i.e. isotropic distribution.

    AFterm, DeltaKQS = geomCalc.deltaLMKQS(EPRXresort, AKQS, phaseConvention = phaseConvention)


    #*** Products
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
    # polProd = (EPRXresort * lambdaTermResort).sum(sumDimsPol)  # Sum polarization terms here to keep total dims minimal in product. Here dims = (mu,mup,Euler/Labels)
    # polProd = (EPRXresort * lambdaTermResort)  # Without polarization terms sum to allow for mupPhase below (reqs. p)
    # Test with alignment term
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



    if SFflagRenorm:
        mTermSumThres.values = mTermSumThres/mTermSumThres.SF

    mTermSumThres['XSraw'] = mTermSumThres.sel({'L':0,'M':0}).drop('LM').copy()  # This basically works, and keeps all non-summed dims... but may give issues later...? Make sure to .copy(), otherwise it's just a pointer.

    # Rescale by sqrt(4pi)*SF, this matches GetCro XS outputs in testing.
    mTermSumThres['XSrescaled'] = mTermSumThres['XSraw']*mTermSumThres['SF']*np.sqrt(4*np.pi)

    mTermSumThres['XSiso'] = XSmatE/3  # ePolyScat defn. for LF cross-section. (i.e. isotropic distribution)
    # mTermSumThres['XS2'] = symDegen * XSmatE/3  # Quick hack for testing, with symDegen

    # Renorm betas by B00?
    if BLMRenorm:
        # mTermSumThres /= mTermSumThres.sel({'L':0,'M':0}).drop('LM')
        if BLMRenorm == 1:
            # Renorm by isotropic XS only
            mTermSumThres /= mTermSumThres['XSraw']

        elif BLMRenorm == 2:
            # Renorm by full t-dependent XS only
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
            mTermSumThres *= normBeta

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
    BetasNormX.attrs['thres'] = thres

    # TODO: update this for **vargs
    # BLMXout.attrs['sumDims'] = sumDims # May want to explicitly propagate symmetries here...?
    # BLMXout.attrs['selDims'] = [(k,v) for k,v in selDims.items()]  # Can't use Xarray to_netcdf with dict set here, at least for netCDF3 defaults.
    BetasNormX.attrs['dataType'] = 'BLM'

    return mTermSumThres, mTermSum, mTerm, BetasNormX
