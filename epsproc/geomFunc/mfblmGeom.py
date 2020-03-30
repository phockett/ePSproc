
import numpy as np

# from epsproc.util import matEleSelector   # Circular/undefined import issue - call in function instead for now.
from epsproc.geomFunc import geomCalc
# from epsproc.geomFunc.geomCalc import (EPR, MFproj, betaTerm, remapllpL, w3jTable,)
from epsproc.geomFunc.geomUtils import genllpMatE

# Code as developed 16/17 March 2020.
# Needs some tidying, and should implement BLM Xarray attrs and format for output.
def mfblmXprod(matEin, QNs = None, EPRX = None, p=[0], BLMtable = None,
                lambdaTerm = None, RX = None, eulerAngs = None,
                thres = 1e-2, thresDims = 'Eke', selDims = {'it':1, 'Type':'L'},
                sumDims = ['mu', 'mup', 'l','lp','m','mp'], sumDimsPol = ['P','R','Rp','p'], symSum = True,
                SFflag = True, squeeze = False, phaseConvention = 'S'):
    """
    Implement :math:`\beta_{LM}^{MF}` calculation as product of tensors.

.. math::
    \begin{eqnarray}
    \beta_{L,-M}^{\mu_{i},\mu_{f}} & = & \sum_{l,m,\mu}\sum_{l',m',\mu'}(-1)^{(\mu'-\mu_{0})}{B_{L,-M}}\nonumber \\
     & \times & \sum_{P,R',R}{E_{P-R}(\hat{e})\Lambda_{R',R}(R_{\hat{n}})}I_{l,m,\mu}^{p_{i}\mu_{i},p_{f}\mu_{f}}(E)I_{l',m',\mu'}^{p_{i}\mu_{i},p_{f}\mu_{f}*}(E)
    \end{eqnarray}


    Where each component is defined by fns. in :py:module:`epsproc.geomFunc.geomCalc` module.

    16/03/20 In progress!

    Parameters
    ----------

    phaseConvention : optional, str, default = 'S'
        Set phase conventions with :py:func:`epsproc.geomCalc.setPhaseConventions`.
        To use preset phase conventions, pass existing dictionary.

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
    if symSum:
        if 'Sym' in matE.dims:
            matE = matE.sum('Sym')  # Sum over ['Cont','Targ','Total'] stacked dims.

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
        if eulerAngs is not None:
            RX = setPolGeoms(eulerAngs = eulerAngs)

        # *** Lambda term
        lambdaTerm, lambdaTable, lambdaD, _ = geomCalc.MFproj(RX = RX, form = 'xarray', phaseConvention = phaseConvention)
        # lambdaTermResort = lambdaTerm.squeeze().drop('l').drop('lp')   # This removes photon (l,lp) dims fully.
        lambdaTermResort = lambdaTerm.squeeze(['l','lp']).drop(['l','lp'])  # Safe squeeze & drop of selected singleton dims only.

    # *** Blm term with specified QNs
    if BLMtable is None:

        QNsBLMtable = QNs.copy()

        # Switch signs (m,M) before 3j calcs.
        if phaseCons['mfblmCons']['BLMmPhase']:
            QNsBLMtable[:,3] *= -1
            QNsBLMtable[:,5] *= -1

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


    #*** Products
    # Matrix element pair-wise multiplication by (l,m,mu) dims
    matEconj = matEthres.copy().conj()
    matEconj = matEconj.unstack().rename({'l':'lp','m':'mp','mu':'mup'})
    matEmult = matEconj * matEthres.unstack()
    matEmult.attrs['dataType'] = 'multTest'

    # Threshold product and drop dims.
    # matEmult = ep.util.matEleSelector(matEmult, thres = thres, dims = thresDims)
    matEmult = matEleSelector(matEmult, thres = thres, dims = thresDims)

    # Product terms with similar dims
    BLMprod = matEmult * BLMtableResort  # Unstacked case with phase correction
    # polProd = (EPRXresort * lambdaTermResort).sum(sumDimsPol)  # Sum polarization terms here to keep total dims minimal in product. Here dims = (mu,mup,Euler/Labels)
    polProd = (EPRXresort * lambdaTermResort)  # Without polarization terms sum to allow for mupPhase below (reqs. p)

    # Set additional phase term, (-1)^(mup-p) **** THIS MIGHT BE SPURIOUS FOR GENERAL EPR TENSOR CASE??? Not sure... but definitely won't work if p summed over above!
    if phaseCons['mfblmCons']['mupPhase']:
        mupPhaseTerm = np.power(-1, np.abs(polProd.mup - polProd.p))
        polProd *= mupPhaseTerm

    polProd = polProd.sum(sumDimsPol)

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

    # Normalise
    # TODO: Set XS as per old mfpad()
#     BLMXout['XS'] = (('Eke','Euler'), BLMXout[0].data)  # Set XS = B00
#     BLMXout = BLMXout/BLMXout.XS  # Normalise
    mTermSumThres['XS'] = mTermSumThres.sel({'L':0,'M':0}).drop('LM').copy()  # This basically works, and keeps all non-summed dims... but may give issues later...? Make sure to .copy(), otherwise it's just a pointer.
    mTermSumThres /= mTermSumThres.sel({'L':0,'M':0}).drop('LM')

    # Propagate attrs
    mTermSum.attrs = mTerm.attrs
    mTermSum.attrs['dataType'] = 'multTest'

    mTermSumThres.attrs = mTerm.attrs
    mTermSumThres.attrs['dataType'] = 'multTest'

    return mTermSumThres, mTermSum, mTerm
