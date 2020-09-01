
import numpy as np

# from epsproc.util import matEleSelector   # Circular/undefined import issue - call in function instead for now.
from epsproc.sphCalc import setPolGeoms
from epsproc.geomFunc import geomCalc, geomUtils
# from epsproc.geomFunc.geomCalc import (EPR, MFproj, betaTerm, remapllpL, w3jTable,)
# from epsproc.geomFunc.geomUtils import genllpMatE

# Code as developed 16/17 March 2020.
# Needs some tidying, and should implement BLM Xarray attrs and format for output.
def lfblmXprod(matEin, QNs = None, EPRX = None, p=[0], BLMtable = None,
                lambdaTerm = None, RX = None, eulerAngs = None,
                thres = 1e-2, thresDims = 'Eke', selDims = {'it':1, 'Type':'L'},
                sumDims = ['mu', 'mup', 'l','lp','m','mp'], sumDimsPol = ['P','R','Rp','p'], symSum = True,
                SFflag = True, squeeze = False, phaseConvention = 'S',
                dlistMatE = ['lp', 'l', 'L', 'mp', 'm', 'M'],
                dlistP = ['p1', 'p2', 'L', 'mup', 'mu', 'M'],
                normFactors = [1/3, 1/5]):
    r"""
    Implement :math:`\beta_{LM}` calculation as product of (Clebsch-Gordan coeff) tensors.

    Formalism as per *Cross section and asymmetry parameter calculation for sulfur 1s photoionization of SF6*, A. P. P. Natalense and R. R. Lucchese, J. Chem. Phys. 111, 5344 (1999), http://dx.doi.org/10.1063/1.479794

    .. math::
        \begin{eqnarray}
        \beta_{\mathbf{k}}^{L,V} & = & \frac{3}{5}\frac{1}{\sum_{p\mu lhv}|I_{\mathbf{k},\hat{n}}^{(L,V)}|^{2}}\sum_{\stackrel{p\mu lhvmm_{v}}{p'\mu'l'h'v'm'm'_{v}}}(-1)^{m'-m_{v}}I_{\mathbf{k},\hat{n}}^{(L,V)} \\
        & \times & \left(I_{\mathbf{k},\hat{n}}^{(L,V)}\right)^{*}b_{lhm}^{p\mu}b_{l'h'm'}^{p'\mu'*}b_{1vm_{v}}^{p_{v}\mu_{v}}b_{1v'm'_{v}}^{p'_{v}\mu'_{v}*} \\
        & \times & [(2l+1)(2l'+1)]^{1/2}(1100|20)(l'l00|20) \\
        & \times & (11-m'_{v}m_{v}|2M')(l'l-m'm|2-M')
        \end{eqnarray}


    Where each component is defined by fns. in :py:module:`epsproc.geomFunc.geomCalc` module.

    22/06/20 In progress - code crudely adapted from mf/af cases, rather messy.
    Dev code:
        http://localhost:8888/lab/tree/dev/ePSproc/geometric_method_dev_Betas_090320.ipynb
        D:\code\ePSproc\python_dev\ePSproc_MFBLM_Numba_dev_tests_120220.PY
        http://localhost:8888/notebooks/epsproc/tests/LF_AF_verification_tests_190620_fullCG.ipynb
    TOTAL MESS AT THE MOMENT>>?>>>>?DFdas<>r ty



    Parameters
    ----------

    phaseConvention : optional, str, default = 'S'
        Set phase conventions with :py:func:`epsproc.geomCalc.setPhaseConventions`.
        To use preset phase conventions, pass existing dictionary.


    normFactors : optional, list of int or float, default = [1/3, 1/5]
        Additional normalization factor to match ePS defns.
        Should be mu and m degeneracy factor? Always 1/3 for B0 term, 1/5 for B2 term?
        Or 1/(2L+1) term?

    """

    from epsproc.util import matEleSelector

    # For transparency/consistency with subfunctions, str/dict now set in setPhaseConventions()
    phaseCons = geomCalc.setPhaseConventions(phaseConvention = phaseConvention)

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
        if 'Sym' in matEthres.dims:
            matEthres = matEthres.sum('Sym')  # Sum over ['Cont','Targ','Total'] stacked dims.


    # Calculate CG sets

    #*** For (l,lp) terms, use code from betaTerm(), but CG version
    # Set QNs
    QNs = geomUtils.genllpMatE(matEthres, phaseConvention = phaseConvention)
    # Set phase convention for term (lp,l,-mp,m|2,-M) - NOTE also ordering of lp,l etc.!

    # dlistMatE = ['lp', 'l', 'L', 'mp', 'm', 'M']  # Now passed.
    if phaseCons['lfblmCGCons']['negmp']:
        QNs[:,3] *= -1  # mp > -mp  # Including this restricts things to mp=0 only? Phase con choice? Doesn't seem correct!
                    # Full calculation including this term sends PU continuum to zero/NaN.

    if phaseCons['lfblmCGCons']['negM']:
        QNs[:,5] *= -1  # M > -M

    # Calculate main term + term with m=mp=M=0
    CGmatE = geomCalc.CG(QNs=QNs, dlist = dlistMatE)
    CGmatEM0 = geomCalc.CG(QNs = geomUtils.genllLList(QNs, uniqueFlag = True, mFlag = False),
                  dlist = dlistMatE, form='xdaLM')

    # Multiply terms (as per betaTerm() code)
    CGmatEmult = CGmatE.unstack()*CGmatEM0.drop('mSet').squeeze().unstack()

    # Reset any flipped coord values for consistency before phase & matE multiplication terms.
    if phaseCons['lfblmCGCons']['negmp']:
        CGmatEmult['mp'] *= -1  # mp > -mp  # Including this restricts things to mp=0 only? Phase con choice? Doesn't seem correct!
                    # Full calculation including this term sends PU continuum to zero/NaN.

    if phaseCons['lfblmCGCons']['negM']:
        CGmatEmult['M'] *= -1  # M > -M


    #*** For photon terms
    # See EPR() and MFproj() for existing 3j prototypes.
    # dlistP = ['p1', 'p2', 'L', 'mup', 'mu', 'M']  # Note p1=p2=1 and will be dropped
    QNsP=geomUtils.genllL(Lmin=1,Lmax=1)  # Generates all possible l1=l2=1 terms, no phase convention applied.

    if phaseCons['lfblmCGCons']['negmup']:
        QNsP[:,3] *= -1  # mup > -mup

    if phaseCons['lfblmCGCons']['negMP']:
        QNsP[:,5] *= -1  # M > -M

    CGP = geomCalc.CG(QNs=QNsP, dlist = dlistP, form='xdaLM')
    CGPM0 = geomCalc.CG(QNs = geomUtils.genllLList(QNsP, uniqueFlag = True, mFlag = False),
                  dlist = dlistP, form='xdaLM')

    # Multiply terms (as per betaTerm() code)
    CGPmult = CGP.unstack()*CGPM0.drop('mSet').squeeze().unstack()

    # Reset any flipped coord values for consistency before phase & matE multiplication terms.
    if phaseCons['lfblmCGCons']['negmup']:
        CGPmult['mup'] *= -1  # mup > -mup

    if phaseCons['lfblmCGCons']['negMP']:
        CGPmult['M'] *= -1  # M > -M

    # *** Multiply all geometric terms
    # Multiply matE * photon terms
    CGPM = CGmatEmult * CGPmult.squeeze(['p1','p2']).drop(['p1','p2'])

    # Set phase = (-1)^(m'-mu)
    CGPM *= np.power(-1, np.abs(CGPM.mp - CGPM.mu))

    #*** Full multiplication by matrix elments.
    # Define matrix element multiplication term (code from mfblmGeom())
    matEconj = matEthres.copy().conj()
    # matEconj = matEconj.unstack().rename({'l':'lp','m':'mp','mu':'mup'})  # Full unstack
    # matEmult = matEthres.unstack() * matEconj
    matEconj = matEconj.unstack('LM').rename({'l':'lp','m':'mp','mu':'mup'})  # Unstack LM only.
    matEmult = matEthres.unstack('LM') * matEconj
    matEmult.attrs['dataType'] = 'multTest'

    # Threshold product and drop dims.
    # matEmult = ep.util.matEleSelector(matEmult, thres = thres, dims = thresDims)
    matEmult = matEleSelector(matEmult, thres = thres, dims = thresDims)

    # Multiplication & sum to BLM values
    degenllp = np.sqrt((2*matEmult.l+1)*(2*matEmult.lp+1))
    Betas = (CGPM * matEmult * degenllp).sum(sumDims).stack({'LM':['L','M']})

    # Additional factors & renorm
    #  XSmatE = (matE * matE.conj()).sel(selDims).sum(['LM','mu'])  # (['LM','mu','it'])  # Cross section as sum over mat E elements squared (diagonal terms only)
    XSmatE = normFactors[0] * (matEthres * matEthres.conj()).sum(['LM','mu'])  # Use selected & thresholded matE.
                                # NOTE - this may fail in current form if dims are missing.
    normBeta = normFactors[1] * (1/XSmatE)  # Normalise by sum over matrix elements squared.

    BetasNorm = Betas.copy()    # Check different renorm methods
    BetasNormXS = normBeta * Betas.copy()

    BetasNorm['XS'] = BetasNorm.sel({'L':0,'M':0}).drop('LM').copy()  # Set XS as B00 term
    # BetasNorm['XS'] = XSmatE  BetasNorm.sel({'L':0,'M':0}).drop('LM').copy()  # Set XS as B00 term
    BetasNorm /= BetasNorm.sel({'L':0,'M':0}).drop('LM')  # Renorm BLM/B00
    BetasNorm *= 1/5  # Renorm by 1/5 == 1/(2L+1)?  Could be included correctly above, or just normalisation choice in ePS?

    BetasNorm['XS2'] = XSmatE

    return BetasNormXS, BetasNorm, Betas, XSmatE
