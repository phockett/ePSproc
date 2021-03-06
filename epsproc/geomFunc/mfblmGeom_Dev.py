# IGNORE THIS - dev version playing only.

def mfblmXprod(matEin, QNs = None, EPRX = None, p=[0], lambdaTerm = None, BLMtable = None,
               thres = 1e-2, thresDims = 'Eke',
               selDims = {'it':1, 'Type':'L'}, sumDims = ['mu', 'mup', 'l','lp','m','mp'], sumDimsPol = ['P','R','Rp','p'], symSum = True,
               SFflag = True):

    # Fudge - set this for now to enforce additonal unstack and phase corrections later.
    BLMtableResort = None

    # Threshold and selection
    # Make explicit copy of data to avoid any overwrite issues
    matE = matEin.copy()
    matE.attrs = matEin.attrs  # May not be necessary with updated Xarray versions

    # Use SF (scale factor)
    # Write to data.values to make sure attribs are maintained. (Not the case for da = da*da.SF)
    if SFflag:
        matE.values = matE * matE.SF

    if symSum:
        matE = matE.sum('Sym')  # Sum over ['Cont','Targ','Total'] stacked dims.

    matEthres = ep.util.matEleSelector(matE, thres = thres, inds = selDims, dims = thresDims)

    # Set terms if not passed to function
    if QNs is None:
        QNs = genllpMatE(matEthres)

    #*** Polarization terms
    if EPRX is None:
        # *** EPR
#         EPRX = geomCalc.EPR(form = 'xarray', p = p).sel({'R-p':0})  # Set for R-p = 0 for p=0 case (redundant coord) - need to fix in e-field mult term!
        # EPRXresort = EPRX.unstack().squeeze().drop('l').drop('lp')  # This removes photon (l,lp) dims fully. Be careful with squeeze() - sends singleton dims to non-dimensional labels.
#         EPRXresort = EPRX.unstack().drop('l').drop('lp')  # This removes photon (l,lp) dims fully, but keeps (p,R) as singleton dims.
#         EPRXresort = EPRX.unstack().squeeze(['l','lp']).drop(['l','lp'])  # Safe squeeze & drop of selected singleton dims only.

#         EPRX = geomCalc.EPR(form = 'xarray', p = p).unstack().sum(['p','R-p'])  # Set for general sum over (p,R-p) terms - STILL need to fix in e-field mult term!
#         EPRX = geomCalc.EPR(form = 'xarray', p = p).unstack().sum('R-p')  # Set for general sum over (p,R-p) terms - STILL need to fix in e-field mult term!
        EPRX = geomCalc.EPR(form = 'xarray', p = p).unstack().sel({'R-p':0}).drop('R-p')
        EPRXresort = EPRX.squeeze(['l','lp']).drop(['l','lp'])  # Safe squeeze & drop of selected singleton dims only.

        Rphase = True
        if Rphase:
            EPRXresort['R'] *= -1

    if lambdaTerm is None:
        # *** Lambda term
        lambdaTerm, lambdaTable, lambdaD, _ = geomCalc.MFproj(form = 'xarray', phaseConvention='E')
#         lambdaTermResort = lambdaTerm.squeeze().drop('l').drop('lp')   # This removes photon (l,lp) dims fully.
        lambdaTermResort = lambdaTerm.squeeze(['l','lp']).drop(['l','lp'])  # Safe squeeze & drop of selected singleton dims only.

    # *** Blm term with specified QNs
    if BLMtable is None:
        BLMmPhase = False
        QNsBLMtable = QNs.copy()
        if BLMmPhase:
            QNsBLMtable[:,3] *= -1
            QNsBLMtable[:,5] *= -1

        BLMtable = geomCalc.betaTerm(QNs = QNsBLMtable, form = 'xdaLM')

#         if BLMmPhase:
#             BLMtable['m'] *= -1

    if BLMtableResort is None:
        # Apply additional phase convention
        BLMtableResort = BLMtable.copy().unstack()
        Mphase = True
        mphase = False

#         if BLMmPhase:
#             BLMtableResort['m'] *= -1  # Sign flip to undo m > -m setting in betaTerm calc.
#             BLMtableResort['M'] *= -1  # Sign flip to undo m > -m setting in betaTerm calc.

        if Mphase:
            BLMtableResort['M'] *= -1
            BLMtableResort *= np.power(-1, np.abs(BLMtableResort.M))  # Associated phase term

        if mphase:
            BLMtableResort['m'] *= -1
            BLMtableResort *= np.power(-1, np.abs(BLMtableResort.m))  # Associated phase term



    #*** Thresholding & dim dropping...
    # To minimize ND size
    # Q: can use matEleSelector here, but might be an issue with large arrays due to copying.
    # Using xr.where() directly should prevent array copy.

    matEmax = matEthres.max()  # Set as a scale factor for later
    thresScaled = matEmax * thres

    BLMtableResort = ep.util.matEleSelector(BLMtableResort, thres=thresScaled)  # Reduce size before big multiplication


    #*** Products
    # Matrix element pair-wise multiplication by (l,m,mu) dims
    matEconj = matEthres.copy().conj()
    matEconj = matEconj.unstack().rename({'l':'lp','m':'mp','mu':'mup'})
    matEmult = matEconj * matEthres.unstack()
    matEmult.attrs['dataType'] = 'multTest'

    # Threshold product and drop dims.
    matEmult = ep.util.matEleSelector(matEmult, thres = thres, dims = thresDims)

    # Product terms with similar dims
    BLMprod = matEmult * BLMtableResort  # Unstacked case with phase correction
    # polProd = (EPRXresort * lambdaTermResort).sum(sumDimsPol)  # Sum polarization terms here to keep total dims minimal in product. Here dims = (mu,mup,Euler/Labels)
    polProd = (EPRXresort * lambdaTermResort)  # Without polarization terms sum to allow for mupPhase below (reqs. p)

    mupPhase = True
    # Set additional phase term, (-1)^(mup-p)  **** THIS MIGHT BE SPURIOUS FOR GENERAL EPR TENSOR CASE??? Not sure... but definitely won't work if p summed over above!
    if mupPhase:
        mupPhaseTerm = np.power(-1, np.abs(polProd.mup - polProd.p))
#         mupPhaseTerm = np.power(-1, np.abs(mTerm.mup - mTerm.p))
#         mupPhaseTerm = np.power(-1, np.abs(mTerm.mup))   # TESTING ONLY, with p already summed over... WITHOUT THIS PU cont (N2 test case) becomes -ve!
        polProd *= mupPhaseTerm

    polProd = polProd.sum(sumDimsPol)


    # Test big mult...
#     mTerm = polProd.sel({'R':0,'Labels':'z'}) * BLMprod.sum(['Total'])    # With selection of z geom.  # BLMprod.sum(['Cont', 'Targ', 'Total'])
#     mTerm = polProd.sel({'R':0}) * BLMprod.sum(['Total'])    # BLMprod.sum(['Cont', 'Targ', 'Total'])
    mTerm = polProd * BLMprod
    # Multiplication works OK, and is fast... but might be an ugly result... INDEED - result large and slow to manipulate, lots of dims and NaNs. Better to sub-select terms first!

    # No subselection, mTerm.size = 6804000
    # For polProd.sel({'R':0}), mTerm.size = 1360800
    # For polProd.sel({'R':0,'Labels':'z'}), mTerm.size = 453600
    # Adding also BLMprod.sum(['Total']), mTerm.size = 226800
    # Adding also BLMprod.sum(['Cont', 'Targ', 'Total']), mTerm.size = 113400  So, for sym specific calcs, may be better to do split-apply type methods

#     mupPhase = True
#     # Set additional phase term, (-1)^(mup-p)  **** THIS MIGHT BE SPURIOUS FOR GENERAL EPR TENSOR CASE??? Not sure... but definitely won't work if p summed over above!
#     if mupPhase:
# #         mupPhaseTerm = np.power(-1, np.abs(mTerm.mup - mTerm.p))
#         mupPhaseTerm = np.power(-1, np.abs(mTerm.mup))   # TESTING ONLY, with p already summed over... WITHOUT THIS PU cont (N2 test case) becomes -ve!
#         mTerm = mTerm * mupPhaseTerm

    mTerm.attrs['file'] = 'MulTest'  # Temporarily adding this, not sure why this is an issue here however (not an issue for other cases...)


    # Sum and threshold
#     sumDims = ['P', 'mu', 'mup', 'Rp', ]  # Define dims to sum over
    xDim = {'LM':['L','M']}
    mTermSum = mTerm.sum(sumDims).squeeze()
    mTermSumThres = matEleSelector(mTermSum.stack(xDim), thres=thres, dims = thresDims)
#     mTermSumThres = mTermSum

    # Normalise
    # TODO: Set XS as per old mfpad()
#     BLMXout['XS'] = (('Eke','Euler'), BLMXout[0].data)  # Set XS = B00
#     BLMXout = BLMXout/BLMXout.XS  # Normalise
    mTermSumThres['XS'] = mTermSumThres.sel({'L':0,'M':0}).drop('LM').copy()  # This basically works, and keeps all non-summed dims. Make sure to .copy(), otherwise it's just a pointer.
    mTermSumThres /=  mTermSumThres['XS']    # mTermSumThres.sel({'L':0,'M':0}).drop('LM')
    mTermSumThres.attrs['dataType'] = 'multTest'

    return mTermSumThres, mTermSum, mTerm
