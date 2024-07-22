"""
ePSproc Gamma functions: calculations

- Set test matrix elements.
- Compute state-resolved betas from legacy Gamma files.

17/07/24

"""

import pandas as pd
import numpy as np

from epsproc.geomFunc.geomCalc import w3jTable


#**** Basic functions for legacy or new gamma calcs.
def setTestMatE(gammaDF, rand=False, phase=True):
    """
    Set test matrix elements to Pandas DF.
    
    Matches (l,m,lam) terms from input gammas.
    
    Default case assigns 0.1*l+0.1
    
    If rand=True, assign random values
    
    If phase=True, assign random phase.
    
    """
    
    matE = []

    # for l1 in gammaDF.index.levels[2]:
        # for lam1 in gammaDF.index.levels[3]:
    
    # Use labels instead of numerical index
    for l1 in gammaDF.index.get_level_values('l1').unique():
        for lam1 in gammaDF.index.get_level_values('lambda1').unique():
            
            if rand:
                llamValue = np.random.rand()
            else:
                llamValue = 0.1+l1*0.1
                
            if phase:
                llamValue = llamValue +1j*np.random.rand()
            
            matE.append([l1,lam1, llamValue])
            

    matEdf = pd.DataFrame(matE, columns=['l1','lambda1','matE1'])
    matEdf.set_index(keys=['l1','lambda1'], inplace=True)
    
    return matEdf
    

def setMatEPrime(matE1):
    """
    Set prime matE from existing DF.
    """
    
    # Set prime terms...
    matEdf2 = matE1.copy()
    matEdf2.index.rename({'l1':'l2','lambda1':'lambda2'},inplace=True)
    matEdf2.rename(columns={'matE1':'matE2'}, inplace=True)
    
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
        
    matE2 = setMatEPrime(matE1)
    
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



def betaCalc(gammaDF, matE=None, **kwargs):
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

def Ccalc(channel=None, lmax=4, thres=1e-4):
    """
    Compute C-params for given channel and lmax.
    
    Method: compute all 3j terms, then subselect & multiply as PD DataFrames.
    
    NOTE - C2/primed terms currently assumes single 'p', also INCOHERENT over Nc/Mc as per legacy codes.
    But may need to revisit this for general case.
    
    Parameters
    ----------
    channel : optional, list or array
        [Ni,Ki,N+,K+]
        
        NOTE: set K=None to run for all allowed terms for given N.
        If not set, run for test case Ni=2, N+=1, all K.
        
    lmax : optional, int, default = 4
    
    
    thres : optional, float or None, default = 1e-4
        Apply threshold to abs(C) product terms, and drop.
        If None, skip thresholding.
        
    """

    # Set master table of 3j results
    pdmaster = w3jTable(Lmax = lmax, form = 'pd', nonzeroFlag = True)
    
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
    tjLFNt = tjLFNt.xs(Ni,level='Ni').xs(Nc,level='Nc')

    # MF Nt coupling
    tjMFNtTerms = {'l':'Nc','lp':'Ni','L':'Nt','m':'Kc','mp':'Ki','M':'Kt'}
    tjMFNt = pdmaster.copy()
    tjMFNt.index.rename(tjMFNtTerms, inplace=True)
    tjMFNt = tjMFNt.xs(Ni,level='Ni').xs(Nc,level='Nc')
    
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
    
    if thres is not None:
        C1 = C1[C1.pipe(np.abs) > thres].dropna()

    #*** Set C2/prime terms
    # NOTE - currently assumes single 'p', also INCOHERENT over Nc/Mc?
    # May need to revisit this for general case.
    C2 = C1.copy()
    C2.rename(columns={'C1':'C2'}, inplace=True)
    C2.index.rename({'l':'lp','m':'mp','lam':'lamp','Nt':'Ntp','Mi':'Mip','q':'qp'}, inplace=True)

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
    
    return Cprod



def gammaCalc(channel=None,Cterms = None, denMat = None, 
              sumList = ['q','qp','Mi','Mip','Nt','Ntp','Mc'],
              **kwargs):
    """
    Compute general gamma parameters.
    
    TODO:
    - Implement BetaTerm (same as general case?)
    - Implement density matrix multiplication
    
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
        
    sumList : list, default = ['q','qp','Mi','Mip','Nt','Ntp','Mc']
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
        print("*** Density matrix mult not yet implemented.")
        
    
    # Sum over some dims
    dims = list(Cterms.index.names)
    # sumList = ['q','qp','Mi','Mip','Nt','Ntp','Mc']

    dimsGroup = list({*dims}-{*sumList})

    gammaPD = Cterms.groupby(by=dimsGroup).sum()

    gammaPD['prod'] = (2*Ni+1) * (2*Nc+1) * \
                      gammaPD['prod']

    gammaPD['lPhase'] = (-1j)**((gammaPD.index.get_level_values(level='lp')-gammaPD.index.get_level_values(level='l')).values)
    gammaPD['prod'] *= gammaPD['lPhase']
    
    gammaPD.rename(columns={'prod':'gamma'}, inplace=True)
    
    return gammaPD