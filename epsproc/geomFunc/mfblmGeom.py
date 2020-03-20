
import numpy as np

from epsproc.util import matEleSelector
from epsproc.geomFunc import geomCalc
# from epsproc.geomFunc.geomCalc import (EPR, MFproj, betaTerm, remapllpL, w3jTable,)
from epsproc.geomFunc.geomUtils import genllpMatE

# Code as developed 16/17 March 2020.
# Needs some tidying, and should implement BLM Xarray attrs and format for output.
def mfblmXprod(matE, QNs = None, EPRX = None, p=[0], lambdaTerm = None, BLMtable = None, thres = 1e-2, thresDims = 'Eke',
               selDims = {'it':1, 'Type':'L'}, sumDims = ['P', 'mu', 'mup', 'Rp','l','lp','m','mp'], squeeze = False):
    """
    Implement :math:`\beta_{LM}^{MF}` calculation as product of tensors.

.. math::
    \begin{eqnarray}
    \beta_{L,-M}^{\mu_{i},\mu_{f}} & = & \sum_{l,m,\mu}\sum_{l',m',\mu'}(-1)^{(\mu'-\mu_{0})}{B_{L,-M}}\nonumber \\
     & \times & \sum_{P,R',R}{E_{P-R}(\hat{e})\Lambda_{R',R}(R_{\hat{n}})}I_{l,m,\mu}^{p_{i}\mu_{i},p_{f}\mu_{f}}(E)I_{l',m',\mu'}^{p_{i}\mu_{i},p_{f}\mu_{f}*}(E)
    \end{eqnarray}


    Where each component is defined by fns. in :py:module:`epsproc.geomFunc.geomCalc` module.

    16/03/20 In progress!

    """


    # Fudge - set this for now to enforce additonal unstack and phase corrections later.
    BLMtableResort = None

    # Threshold and selection
    matEthres = matEleSelector(matE, thres = thres, inds = selDims, dims = thresDims)

    # Set terms if not passed to function
    if QNs is None:
        QNs = genllpMatE(matEthres)

    #*** Polarization terms
    if EPRX is None:
        # *** EPR
        EPRX = geomCalc.EPR(form = 'xarray', p = p).sel({'R-p':0})  # Set for R-p = 0 for p=0 case (redundant coord) - need to fix in e-field mult term!
        # EPRXresort = EPRX.unstack().squeeze().drop('l').drop('lp')  # This removes photon (l,lp) dims fully. Be careful with squeeze() - sends singleton dims to non-dimensional labels.
        # EPRXresort = EPRX.unstack().drop('l').drop('lp')  # This removes photon (l,lp) dims fully, but keeps (p,R) as singleton dims.
        EPRXresort = EPRX.unstack().squeeze(['l','lp']).drop(['l','lp'])  # Safe squeeze & drop of selected singleton dims only.

        Rphase = True
        if Rphase:
            EPRXresort['R'] *= -1

    if lambdaTerm is None:
        # *** Lambda term
        lambdaTerm, lambdaTable, lambdaD, _ = geomCalc.MFproj(form = 'xarray')
        # lambdaTermResort = lambdaTerm.squeeze().drop('l').drop('lp')   # This removes photon (l,lp) dims fully.
        lambdaTermResort = lambdaTerm.squeeze(['l','lp']).drop(['l','lp'])  # Safe squeeze & drop of selected singleton dims only.

    # *** Blm term with specified QNs
    if BLMtable is None:
        BLMtable = geomCalc.betaTerm(QNs = QNs, form = 'xdaLM')

    if BLMtableResort is None:
        # Apply additional phase convention
        BLMtableResort = BLMtable.copy().unstack()
        Mphase = True
        mphase = True

        if Mphase:
            BLMtableResort['M'] *= -1
            BLMtableResort *= np.power(-1, np.abs(BLMtableResort.M))  # Associated phase term

        if mphase:
            BLMtableResort['m'] *= -1
            BLMtableResort *= np.power(-1, np.abs(BLMtableResort.m))  # Associated phase term


    #*** Products
    # Matrix element pair-wise multiplication by (l,m,mu) dims
    matEconj = matEthres.copy().conj()
    matEconj = matEconj.unstack().rename({'l':'lp','m':'mp','mu':'mup'})
    matEmult = matEconj * matEthres.unstack()
    matEmult.attrs['dataType'] = 'multTest'

    # Product terms with similar dims
    BLMprod = matEmult * BLMtableResort  # Unstacked case with phase correction
    polProd = EPRXresort * lambdaTermResort

    # Test big mult...
    # mTerm = polProd.sel({'R':0,'Labels':'z'}) * BLMprod.sum(['Total'])    # With selection of z geom.  # BLMprod.sum(['Cont', 'Targ', 'Total'])
    mTerm = polProd.sel({'R':0}) * BLMprod    # BLMprod.sum(['Cont', 'Targ', 'Total'])
    # Multiplication works OK, and is fast... but might be an ugly result... INDEED - result large and slow to manipulate, lots of dims and NaNs. Better to sub-select terms first!

    # No subselection, mTerm.size = 6804000
    # For polProd.sel({'R':0}), mTerm.size = 1360800
    # For polProd.sel({'R':0,'Labels':'z'}), mTerm.size = 453600
    # Adding also BLMprod.sum(['Total']), mTerm.size = 226800
    # Adding also BLMprod.sum(['Cont', 'Targ', 'Total']), mTerm.size = 113400  So, for sym specific calcs, may be better to do split-apply type methods

    mupPhase = True
    # Set additional phase term, (-1)^(mup-p)
    if mupPhase:
        mupPhaseTerm = np.power(-1, np.abs(mTerm.mup - mTerm.p))
        mTerm = mTerm * mupPhaseTerm

    # mTerm.attrs['file'] = 'MulTest'  # Temporarily adding this, not sure why this is an issue here however (not an issue for other cases...)
    mTerm.attrs = matE.attrs  # Propagate attrs from input matrix elements.


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
    mTermSumThres.attrs['dataType'] = 'multTest'

    return mTermSumThres, mTermSum, mTerm
