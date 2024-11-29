"""
ePSproc Gamma functions: calculations

- Set test matrix elements.
- Compute state-resolved betas from legacy Gamma files.
- Compute C values (python version)
- Compute state-resolved betas (python version).

17/07/24

NOTE: need to check PD NaN propagation. pd.prod(axis=1) may ignore NaN?  But other pd.sum() methods propagate NaNs - need to check carefully, maybe replace with 0 in cases where there are issues.

"""

import pandas as pd
import numpy as np
import xarray as xr

from epsproc import multiDimXrToPD
from epsproc.geomFunc.geomCalc import w3jTable
from epsproc.geomFunc.geomCalc import betaTerm as betaTermCalc


#**** Basic functions for legacy or new gamma calcs.
def setTestMatE(gammaDF, rand=False, phase=True, cols=['l1','lambda1','matE1'], colsPrime=None):
    """
    Set test matrix elements to Pandas DF.
    
    Matches (l,m,lam) terms from input gammas. (Set `cols=['l1','lambda1','matE1']` to change naming schema.)
    
    Default case assigns 0.1*l+0.1
    
    If rand=True, assign random values
    
    If phase=True, assign random phase.
    
    """
    
    matE = []

    # for l1 in gammaDF.index.levels[2]:
        # for lam1 in gammaDF.index.levels[3]:
    
    # Use labels instead of numerical index
    for l1 in gammaDF.index.get_level_values(cols[0]).unique():
        for lam1 in gammaDF.index.get_level_values(cols[1]).unique():
            
            if rand:
                llamValue = np.random.rand()
            else:
                llamValue = 0.1+l1*0.1
                
            if phase:
                llamValue = llamValue +1j*np.random.rand()
            
            matE.append([l1,lam1, llamValue])
            

    matEdf = pd.DataFrame(matE, columns=cols)
    matEdf.set_index(keys=cols[0:2], inplace=True)
    
    return matEdf
    

def setMatEPrime(matE1, cols=['l1','lambda1','matE1'], colsPrime=['l2','lambda2','matE2']):
    """
    Set prime matE from existing DF.
    """
    
    # Set prime terms...
    matEdf2 = matE1.copy()
    matEdf2.index.rename({cols[0]:colsPrime[0],cols[1]:colsPrime[1]},inplace=True)
    matEdf2.rename(columns={cols[2]:colsPrime[2]}, inplace=True)
    
    return matEdf2
    


def assignMatE(gammaDF, matE=None, **kwargs):
    """
    Assign matrix elements as columns in Pandas DataFrame of gamma values.
    
    If matE=None, use :py:func:`setTestMatE`, and passed **kwargs.
    
    """
    
    if matE is None:
        matE1 = setTestMatE(gammaDF, **kwargs)
    else:
        matE1=matE
        
    matE2 = setMatEPrime(matE1, **kwargs)
    
    # Assign terms via merge
    # Brief: https://stackoverflow.com/a/55366715
    # Docs: https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.merge.html
    df1 = gammaDF.merge(matE1, 
             left_index=True, 
             right_index=True, 
             how='left')

    df2 = df1.merge(matE2, 
             left_index=True, 
             right_index=True, 
             how='left')

    return df2, matE1, matE2


def betaCalc(gammaDF, matE=None, betaTerm=None, returnType = 'beta', dropnan = True, **kwargs):
    """
    Compute betas from gamma values & matrix elements for state-resolved case.
    
    For legacy case pass args to betaCalcLegacy().
    
    TODO:
    - Better dim name handling (currently hard-coded).
    
    Formalism:
    
    $$
    \begin{eqnarray}
    \beta_{L,M}(k) & = & \sum_{ll'}\sum_{\lambda\lambda'}\sum_{mm'}(-1)^{m}\sqrt{\frac{(2l+1)(2l'+1)(2L+1)}{4\pi}}\nonumber \\
     & \mathsf{x} & \left(\begin{array}{ccc}
    l & l' & L\\
    m & -m' & M
    \end{array}\right)\left(\begin{array}{ccc}
    l & l' & L\\
    0 & 0 & 0
    \end{array}\right)\gamma_{\alpha\alpha_{+}l\lambda ml'\lambda'm'}\nonumber \\
     & \mathsf{x} & \boldsymbol{r}_{kl\lambda}\boldsymbol{r}_{kl'\lambda'}e^{i(\eta_{l\lambda}(k)-\eta_{l'\lambda'}(k))}\label{eq:beta-gamma-general}
    \end{eqnarray}
    $$
    
    Here use:
    - gamma values from this code (legacy or python).
        - Legacy files include betaTerm() values.
        - Python version per `gammaCalc()`.
        - Python version uses `geomCalc.betaTerm()` for additional terms above.
    - Matrix elements as passed, or assigned as per options to `assignMatE()`.
    
    23/08/24 added dropna option here. For cases with NaN matE, skipping this may lead to all-NaN outputs.

    """

    #*** For legacy case (from file), use old function
    if ('legacyGamma' in gammaDF.attrs.keys()) and gammaDF.attrs['legacyGamma']:
        dfSum, dfCalc = betaCalcLegacy(gammaDF, matE=None, **kwargs)
        return dfSum, dfCalc 
    
    
    #*** For new case (python calcs)
    # This should handle arb cols in input DF
    
    # Assign matE
    dfCalc, matE1, matE2 = assignMatE(gammaDF,matE, cols=['l','m','matE1'], colsPrime=['lp','mp','matE2'])
    dfMult = dfCalc.iloc[:, :-2].mul(dfCalc['matE1'],axis=0).mul(dfCalc['matE2'],axis=0)
    
    # Compute betaTerm if not precalculated
    if betaTerm is None:
        Lmax = gammaDF.index.get_level_values(level='l').max()
        BLMtable = betaTermCalc(Lmax = Lmax, form = 'xdaLM') 
        
        # Push to DF with new axis - OK if squeeze=False set!
        # This produces single col, multindex output
        betaTerm,_ = multiDimXrToPD(BLMtable.expand_dims(['BLM']), colDims='BLM', squeeze=False)
        
        betaTerm.rename(columns={0:'betaTerm'}, inplace=True)
    
    
    # Mult by betaTerm
    # 23/08/24 added dropna option here. For cases with NaN matE, skipping this may lead to all-NaN outputs.
    #          TODO: more checks/tests here, haven't carefully verified drop is OK.
    if dropnan:
        dfMult = dfMult.dropna(how='all')
        
    # Merge and multiply as per previous cases...
    BLMprod = dfMult.merge(betaTerm, 
             left_index=True, 
             right_index=True,
             how='left')       # Use index from gamma/dfMult as primary.

    # Multiply & sum
    sumTerms = BLMprod.iloc[:, :-1].mul(BLMprod['betaTerm'],axis=0) 
    betaOut = sumTerms.groupby(by=['L','M']).sum()
    betaOutNorm = betaOut/betaOut.loc[0,0]
    
    if returnType == "full":
        return locals()
    
    else:
        return betaOut, betaOutNorm
    
    
    
def betaCalcLegacy(gammaDF, matE=None, **kwargs):
    """
    Compute betas from (legacy) gamma terms and matrix elements.
    
    All values must be as Pandas DataFrames
    
    If matE=None, use :py:func:`setTestMatE`, and passed **kwargs. 
    
    Note initial testing with outputs from NH3 gamma code `ion_rot_gamma_nh3_4d_NS.c`, which includes 1-photon density matrix and symmetry selection rules in output.
    
    For other codes modifications may be required, e.g. including density matrix.
    
    """
    
    # Assign all terms to master DF
    dfCalc, matE1, matE2 = assignMatE(gammaDF,matE, **kwargs)
    
    # Compute product terms
    # Multiply cols
    dfCalc['BLMprod'] = dfCalc['betaTerm']*dfCalc['gamma']*dfCalc['matE1']*dfCalc['matE2']
    
    # Sum terms
    dfSum = dfCalc.groupby(by=['L','M']).sum()
    
    # Norm vals
    dfSum['BLMnorm'] = dfSum['BLMprod']/dfSum['BLMprod'].loc[0,0]
    
    # Propagate attrs
    dfCalc.attrs = gammaDF.attrs.copy()
    dfCalc.attrs['matE1']=matE1
    dfCalc.attrs['matE2']=matE2

    dfSum.attrs = gammaDF.attrs.copy()
    dfSum.attrs['matE1']=matE1
    dfSum.attrs['matE2']=matE2
    
    
    return dfSum, dfCalc



#***** Compute gammas (new style)

def Ccalc(channel=None, lmax=4, halfIntFlag = False, thres=1e-4, spinWeightings = None):
    """
    Compute C-params for given channel and lmax.
    
    Method: compute all 3j terms, then subselect & multiply as PD DataFrames.
    
    NOTE - C2/primed terms currently assumes single 'p', also INCOHERENT over Nc/Mc as per legacy codes.
    But may need to revisit this for general case.
    See also note on channels and Kt below - may need to add coherence here too.
    
    NOTE - spin terms currently NOT implemented. (I.e. J+ == N+)
    UPDATE Nov. 2024: spin terms now implemented via :py:func:`spinWeightings()`, pass result to Ccalc to incorporate spin.
    
    Formalism:
    
    $$
    \begin{eqnarray}
    C(lm\lambda N_{t}M_{i}\mu_{\lambda}) & = & (2N_{t}+1)(-1)^{M_{+}+q}\left(\begin{array}{ccc}
    N_{t} & 1 & l\\
    M_{t} & p & m
    \end{array}\right)\left(\begin{array}{ccc}
    N_{+} & N_{i} & N_{t}\\
    -M_{+} & M_{i} & M_{t}
    \end{array}\right)\nonumber \\
     & \mathsf{x} & \left(\begin{array}{ccc}
    N_{+} & N_{i} & N_{t}\\
    -K_{+} & K_{i} & K_{t}
    \end{array}\right)\left(\begin{array}{ccc}
    N_{t} & 1 & l\\
    -K_{t} & q & -\lambda
    \end{array}\right)\nonumber \\
     & \mathsf{x} & \left(\begin{array}{ccc}
    N_{+} & J_{+} & S_{+}\\
    M_{+} & M_{J+} & M_{S+}
    \end{array}\right)\left(\begin{array}{ccc}
    N_{+} & J_{+} & S_{+}\\
    K_{+} & P_{+} & \Sigma_{+}
    \end{array}\right)\label{eq:geom-params-C}
    \end{eqnarray}
    $$
    
    
    Parameters
    ----------
    channel : optional, list or array
        [Ni,Ki,N+,K+]
        
        NOTE: set K=None to run for all allowed terms for given N.
        If not set, run for test case Ni=2, N+=1, all K.
        Note Kt is currently NOT coherently summed, Kt=Ki-K+ only for state-resolved cases.
        
        To run all cases to lmax, set channel = 4*[None]
        Note this might produce large output.
        
    lmax : optional, int, default = 4
    
    halfIntFlag : bool, optional, default = False
        If True, include 1/2-int terms in QN creation routine.
    
    thres : optional, float or None, default = 1e-4
        Apply threshold to abs(C) product terms, and drop.
        If None, skip thresholding.
    
    spinWeightings : optional, default = None
        Pass outputs from :py:func:`spinWeightings()` to include spin weightings.
    
    
    Notes
    -----
    27/11/24: added additional xs handling to allow for all Ni,Nc case
              added spinWeightings option, note this currently needs to be run separately and passed as a Pandas DataFrame, as output by :py:func:`spinWeightings()`
    
    07/08/24: tidying up. 
              - Added missing phase and degen factors.
              - Updated docs
              
    NOTE: need to check PD NaN propagation. pd.prod(axis=1) may ignore NaN?  But other pd.sum() methods propagate NaNs - need to check carefully, maybe replace with 0 in cases where there are issues.
    
    """

    # Set master table of 3j results
    pdmaster = w3jTable(Lmax = lmax, form = 'pd', nonzeroFlag = True, halfIntFlag = halfIntFlag)
    
    # Test case, run for all K.
    if channel is None:
        channel = [2,None,1,None]
    
    Ni,Ki,Nc,Kc=channel
    
    #**** Calc all tj terms, then multiply
    #
    # TODO: work out correct merge/join/multiply/cross-prod here... not sure if missing terms currently.
    # UPDATE: sensible way is to index by tj with coupled terms defined/restricted, so always a subset?
    
    # LF photon coupling
    tjLFlTerms = {'l':'Nt','lp':'P','L':'l','m':'Mt','mp':'p','M':'m'}
    tjLFl = pdmaster.copy()
    tjLFl.index.rename(tjLFlTerms, inplace=True)
    tjLFl = tjLFl.xs(1,level='P')  #.xs(Nc,level='Nc')

    # MF photon coupling
    tjMFlTerms = {'l':'Nt','lp':'P','L':'l','m':'Kt','mp':'q','M':'lam'}
    tjMFl = pdmaster.copy()
    tjMFl.index.rename(tjMFlTerms, inplace=True)
    tjMFl = tjMFl.xs(1,level='P')  #.xs(Nc,level='Nc')

    # LF Nt coupling
    tjLFNtTerms = {'l':'Nc','lp':'Ni','L':'Nt','m':'Mc','mp':'Mi','M':'Mt'}
    tjLFNt = pdmaster.copy()
    tjLFNt.index.rename(tjLFNtTerms, inplace=True)
    # tjLFNt = tjLFNt.xs(Ni,level='Ni').xs(Nc,level='Nc')

    # MF Nt coupling
    tjMFNtTerms = {'l':'Nc','lp':'Ni','L':'Nt','m':'Kc','mp':'Ki','M':'Kt'}
    tjMFNt = pdmaster.copy()
    tjMFNt.index.rename(tjMFNtTerms, inplace=True)
    # tjMFNt = tjMFNt.xs(Ni,level='Ni').xs(Nc,level='Nc')
    
    # Subselect on Ni, Nc if set
    if Ni is not None:
        tjLFNt = tjLFNt.xs(Ni,level='Ni')
        tjMFNt = tjMFNt.xs(Ni,level='Ni')
    if Nc is not None:
        tjLFNt = tjLFNt.xs(Nc,level='Nc')
        tjMFNt = tjMFNt.xs(Nc,level='Nc')
    
    # Subselect on Ki,Kc if set
    if Ki is not None:
        tjMFNt = tjMFNt.xs(Ki,level='Ki')
    if Kc is not None:
        tjMFNt = tjMFNt.xs(Kc,level='Kc')
        
    
    # Assign terms for products.
    # This should span all allowed subset of terms
    # If NaNs appear then it may indicate issues with indexing/assignments above
    dfprodLF = tjLFNt.merge(tjLFl, 
             left_index=True, 
             right_index=True,
             how='left')

    dfprodMF = tjMFNt.merge(tjMFl, 
             left_index=True, 
             right_index=True,
             how='left')  

    # Product terms
    dfprod = dfprodLF.merge(dfprodMF, 
         left_index=True, 
         right_index=True,
         # how='right')  # 8080 rows, no Nans
         how='left')  # Same result OK

    dfprod['prod']= dfprod.prod(axis=1)

    # Also set LF, MF products for testing
    dfprodLF['prod']= dfprodLF.prod(axis=1)
    dfprodMF['prod']= dfprodMF.prod(axis=1)
    
    # Set C-terms and multiply
    # thres=1e-2
    C1 = dfprod['prod'].to_frame()
    C1.rename(columns={'prod':'C1'}, inplace=True)
    
    # Degen factor 2Nt+1
    degen = (2*C1.index.get_level_values(level='Nt').values) + 1
    C1 =  C1.multiply(degen, axis=0)  # Need to force row-wise multiply here in general.
    
    # Phase factors (-1)^(Mc+q) - note force to positive powers only
    MpqPhase = (-1)**np.abs((C1.index.get_level_values(level='Mc')+C1.index.get_level_values(level='q')).values)
    C1 =  C1.multiply(MpqPhase, axis=0)  # Need to force row-wise multiply here in general.
    
    if thres is not None:
        C1 = C1[C1.pipe(np.abs) > thres].dropna()

    #*** Set C2/prime terms
    # NOTE - currently assumes single 'p', also INCOHERENT over Nc/Mc?
    # May need to revisit this for general case.
    C2 = C1.copy()
    C2.rename(columns={'C1':'C2'}, inplace=True)
    C2.index.rename({'l':'lp','m':'mp','lam':'lamp','Nt':'Ntp','Mt':'Mtp','Mi':'Mip','q':'qp'}, inplace=True)

    Cprod = C2.merge(C1, 
             left_index=True, 
             right_index=True,
             how='right')  
             # how='left') 

    # Quick tests for Ni=2, Nc=1 (all K)
    # thres = 1e-4 >>> 400508 rows
    # thres = 1e-2 >>> 353 rows
    # No NaNs.
    # No change with left or right index merge.
    # Looks OK.
    
    Cprod['prod'] = Cprod.prod(axis=1)
    
    if thres is not None:
        Cprod = Cprod[Cprod['prod'].pipe(np.abs) > thres].dropna()
    
    # If spin weightings are passed, run additional product and clean-up steps.
    if spinWeightings is not None:
        
        # If passed as dict, use 'sum' term
        if isinstance(spinWeightings, dict):
            spinW = spinWeightings['sum']
        else:
            spinW = spinWeightings
            
        # Assign to column
        spinW.rename(columns={'prod':'prodSpin'}, inplace=True)   #.xs(J, level='Nc').rename(columns={'prod':'prodSpin'})
        dfprodSpin = Cprod.merge(spinW, 
             left_index=True, 
             right_index=True,
             # how='left')      
             how='right')
        
        # Spin weighted product
        dfprodSpin['spinWeighted'] = dfprodSpin['prod']*dfprodSpin['prodSpin']
    
        # Rename coloumns to use 'prod' (maybe assumed in later functions)
        dfprodSpin.rename(columns={'prod':'prodUnweighted'}, inplace=True) 
        dfprodSpin.rename(columns={'spinWeighted':'prod'}, inplace=True)
    
        if thres is not None:
            dfprodSpin = dfprodSpin[dfprodSpin['prod'].pipe(np.abs) > thres].dropna()
        
        return dfprodSpin
        
    else:
        return Cprod



def gammaCalc(channel=None,Cterms = None, denMat = None, 
              sumList = ['q','qp','Mi','Mip','Nt','Ntp','Mt','Mtp','Mc'],
              **kwargs):
    """
    Compute general gamma parameters.
    
    TODO:
    - Implement BetaTerm (same as general case?)
    - Implement density matrix multiplication
        - 25/07/24: implemented pmm multiplication. Note that this currently runs `denMatReformat` for Xarray inputs, and doesn't use `channel` specs currently.
    - 29/11/24: updated to handle spin weighted version of Cterms, although currently needs some manual effort.
        
    Formalism:
    
    $$
    \begin{eqnarray}
    \gamma_{\alpha\alpha_{+}l\lambda ml'\lambda'm'} & = & (2N_{i}+1)(2N_{+}+1)(-i)^{l'-l}\sum_{M_{+}}\sum_{M_{i}M_{i}'}\sum_{N_{t}N_{t}'}\sum_{\mu_{\lambda}\mu_{\lambda}'}{}^{J_{i}K_{i}}\boldsymbol{\rho}_{M_{i}M_{i}'}\nonumber \\
     & \mathsf{x} & C(lm\lambda N_{t}M_{i}q)C(l'm'\lambda'N_{t}'M_{i}'q')\label{eq:gamma-state}
    \end{eqnarray}
    $$
    
    Parameters
    ----------

    channel : optional, list or array, default = None
        [Ni,Ki,N+,K+]
        
        NOTE: set K=None to run for all allowed terms for given N.
        If not set, run for test case Ni=2, N+=1, all K.
        
    Cterms : optional, PD DataFrame
        C values to use, as output by Ccalc()
        If not set, Ccalc(**kwargs) will be run.
        
    denMat : optional, PD DataFrame
        Initial state density matrix p(mi,mi').
        If not set, assume all terms are = 1.
        
    sumList : list, default = ['q','qp','Mi','Mip','Nt','Ntp','Mt','Mtp','Mc']
        QNs to sum over in output.
        Default case matches legacy codes.
        Assumes single (Ni,Nc) state, and selected 'p'.
        
    **kwargs
        Additional args passed to Ccalc() if run.
    
    """
    
    # Test case, run for all K.
    # TODO: repeats Ccalc code...
    if channel is None:
        channel = [2,None,1,None]
    
    Ni,Ki,Nc,Kc=channel
    
    if Cterms is None:
        Cterms = Ccalc(channel,**kwargs)     
        
    if denMat is not None:
        # print("*** Density matrix mult not yet implemented.")
        if isinstance(denMat, xr.core.dataarray.DataArray):
            denMat = denMatReformat(denMat)
          
        # Assign and multiply density matrix & gamma terms
        CtermsRho = Cterms.merge(denMat, 
                 left_index=True, 
                 right_index=True,
                 how='left')
        
        # Crho['prod']*=Crho['rho']  # For single col case

        # For multi-col case, plus return clean DF
        # Per https://stackoverflow.com/a/46779778
        # 29/11/24: added cases for with/without spin (additional Cterms in with spin case).
        if 'prodSpin' in Cterms.columns:
            Cpmm = CtermsRho.iloc[:, 7:].mul(CtermsRho['prod'],axis=0)
        else:
            Cpmm = CtermsRho.iloc[:, 3:].mul(CtermsRho['prod'],axis=0)  #.combine_first(Crho)  # Add this to return original vals too.
    
    # For null density matrix, just set product term - this saves additional axis checks etc....?
    # OR may want to keep C1,C1...?
    else:
        Cpmm = pd.DataFrame(Cterms['prod'])
        # Cpmm = Cterms['prod']
    

    # 28/11/24: modified to allow for multiple Ni,Nc terms
    # TODO: move to pre-sum, and use terms in Pandas index if set
    # 29/11/24: moved and set to use index if required.
    # TODO: seems to be working, but may loose column names here?
    if Ni is not None:
        Cpmm *= (2*Ni+1)
    else:
        Cpmm = Cpmm.multiply(2*Cpmm.index.get_level_values(level='Ni')+1, axis = 0)

    if Nc is not None:
        Cpmm *= (2*Nc+1)
    else:
        Cpmm = Cpmm.multiply(2*Cpmm.index.get_level_values(level='Nc')+1, axis = 0)
    
    # Sum over some dims
    gammaPD = sumPDGroups(Cpmm, sumDims=sumList)

    lPhase = (-1j)**((gammaPD.index.get_level_values(level='lp')-gammaPD.index.get_level_values(level='l')).values)
    
    # Multiply by degen & phase factors
    gammaPD =  gammaPD.multiply(lPhase, axis=0)  # Need to force row-wise multiply here in general.

    
    # TODO: Renorm...?
    
    
    # Initial test version - assumes single col gamma only.
#     # Sum over some dims
#     dims = list(Cterms.index.names)
#     # sumList = ['q','qp','Mi','Mip','Nt','Ntp','Mc']

#     dimsGroup = list({*dims}-{*sumList})

#     gammaPD = Cterms.groupby(by=dimsGroup).sum()

#     gammaPD['prod'] = (2*Ni+1) * (2*Nc+1) * \
#                       gammaPD['prod']

#     gammaPD['lPhase'] = (-1j)**((gammaPD.index.get_level_values(level='lp')-gammaPD.index.get_level_values(level='l')).values)
#     gammaPD['prod'] *= gammaPD['lPhase']
    
#     gammaPD.rename(columns={'prod':'gamma'}, inplace=True)
        
    
    # Force to DataFrame for single t case
    if isinstance(gammaPD, pd.core.series.Series):
        gammaPD = pd.DataFrame(gammaPD)
    
    # Set metadata
    gammaPD.attrs['dataType']="gamma"
    gammaPD.attrs['denMat']=True if denMat is not None else False
    gammaPD.attrs['source']="gammaCalc.gammaCalc()"
    gammaPD.attrs['legacyGamma']=False  # 09/08/24 Set this to allow quick switch on legacy gamma from file vs. new python calcs in ancillary functions.
    
    
    return gammaPD, lPhase, Cpmm, Cterms


def denMatReformat(denMat, dimMap = {'M':'Mi','Mp':'Mip'}, 
                   colDims = 't', sumDims = 'default'):
    """
    Reformat density matrix Xarray to PD version for gammaCalc().
    
    NOTE: default config assumes density matrix format from :py:func:`ep.calc.density.densityFromSphTensor()`.
    
    For other cases manual reformatting may be required.
    
    """
    
    if sumDims is not None:
        if sumDims == 'default':
            denMat = denMat.sum(['K','Q'])
        else:
            denMat = denMat.sum(sumDims)
    
    denMatPD, denMatRestack = multiDimXrToPD(denMat, colDims=colDims)
    
    if dimMap is not None:
        denMatPD.index.rename(dimMap,inplace=True)
        
    # TODO: add cleanup and other selectors?
    # # Remove spurious M states?
    # from epsproc.sphFuncs.sphConv import cleanLMcoords, checkSphDims
    # pmmClean = cleanLMcoords(pmm, refDims=['J','M'])
    # pmmClean = cleanLMcoords(pmmClean, refDims=['Jp','Mp'])
        
    return denMatPD


def sumPDGroups(dataIn, sumDims = None):
    """
    Group by & sum PD DataFrame by sumDims only.
    
    """
    
    dims = list(dataIn.index.names)
    
    # Force list for singleton case, otherwise dimsGroup will not be set correctly.
    if not isinstance(sumDims, list):
        sumDims = [sumDims]

    dimsGroup = list({*dims}-{*sumDims})

    return dataIn.groupby(by=dimsGroup).sum()
    
    
    
def spinWeightings(lmax = 3, Sc = 0.5,
                  sumTerms = ['sigSc', 'Msc', 'Mjc', 'Pc'],
                  selectors = None, query = None):
    """
    Compute spin-coupling terms for J states, as used in photoionization calculations:
    
    $$
    \begin{eqnarray}
    \left(\begin{array}{ccc}
    N_{+} & J_{+} & S_{+}\\
    M_{+} & M_{J+} & M_{S+}
    \end{array}\right)\left(\begin{array}{ccc}
    N_{+} & J_{+} & S_{+}\\
    K_{+} & P_{+} & \Sigma_{+}
    \end{array}\right)\label{eq:geom-params-C}
    \end{eqnarray}
    $$
    
    Where
    
    - S = 1/2
    - N = integer.
    - J = 1/2-int terms inc. spin.
    
    Note params are labelled by prefix 'c'(ore) in output, e.g. $N_{+} = Nc$ etc.
    Labelling matches main gamma calculation.
    
    Parameters
    ----------
    
    lmax : int, optional, default=3
        Max ang. mom. to use for tabulations.
        
    Sc : float, optional, default = 0.5
        Default case set for spin decoupling, but can be set to other values if required.
        
    sumTerms : list, optional, default = ['sigSc', 'Msc', 'Mjc', 'Pc']
        Terms to sum over.
        
    selectors : dict, optional, default = None
        Dictionary to subselect terms using pd.xs().
        E.g. {'Nc':1} to select Nc=1 terms.
        Note .xs() should support multiple states as tuples, e.g. {'Nc':(1,2)}, but this fails in testing in PD v1.5.3
        Alternatively, skip selectors here and use mapping on returned dataframe, 
        e.g. dfprodSpin.loc[map(lambda x: x in [1,2], dfprodSpin.index.get_level_values('Nc'))]
        (solution from https://stackoverflow.com/a/77176420)
        
    query : list, optional, default = None
        List of strings to use for pd.query().
        E.g. ['Nc%2==0'] to select even Nc terms only.
    
    
    Returns
    -------
    
    dict
        Contains Pandas tabulations of results.
        - 'full' complete tabulation.
        - 'sub' subselected terms.
        - 'sum' summed terms (from subselection).
    
    """
    
    # Calculate 3j terms
    # lmax = 3
    pdmasterSpin = w3jTable(Lmax = lmax, form = 'pd', nonzeroFlag = True, halfIntFlag=True)
    # pdmasterSpin

    # LF spin coupling
    tjLFspinTerms = {'l':'Nc','lp':'Jc','L':'Sc','m':'Mc','mp':'Mjc','M':'Msc'}
    tjLFspin = pdmasterSpin.copy()
    tjLFspin.index.rename(tjLFspinTerms, inplace=True)
    tjLFspin = tjLFspin.query('Nc%1 == 0').xs(Sc,level='Sc') #.xs(0.5,level='Sc') #.xs(Jc,level='Jc')

    # MF spin coupling
    tjMFspinTerms = {'l':'Nc','lp':'Jc','L':'Sc','m':'Kc','mp':'Pc','M':'sigSc'}
    tjMFspin = pdmasterSpin.copy()
    tjMFspin.index.rename(tjMFspinTerms, inplace=True)
    tjMFspin = tjMFspin.query('Nc%1 == 0').xs(Sc,level='Sc') #.xs(0.5,level='Sc') #.xs(Jc,level='Jc')  #.xs(Nc,level='Nc')
    
    # Assign terms for products.
    # This should span all allowed subset of terms
    # If NaNs appear then it may indicate issues with indexing/assignments above
    dfprodSpin = tjLFspin.merge(tjMFspin, 
         left_index=True, 
         right_index=True,
         # how='left')      
         how='right')
    
    dfprodSpin.rename(columns={'3j_x':'LFspin','3j_y':'MFspin'}, inplace=True)
    
    # Add product term
    # dfprodSpin['prod']=dfprodSpin['3j_x']*dfprodSpin['3j_y']
    dfprodSpin['prod']=dfprodSpin.prod(axis=1)
    
    # Subselect if passed
    # if selectors is None:
    dfprodSub = dfprodSpin
    if selectors is not None:
        for k,v in selectors.items():
            dfprodSub = dfprodSub.xs(v, level=k)
        
    
    # Sum terms - use existing wrapper for this (need to group then sum)
    # For sum by multindex group
    dfprodSum = sumPDGroups(dfprodSub,sumTerms)
    
    return {'sum':dfprodSum, 'sub':dfprodSub, 'full':dfprodSpin}
