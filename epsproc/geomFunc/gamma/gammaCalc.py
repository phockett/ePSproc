"""
ePSproc Gamma functions: IO

- Set test matrix elements.
- Compute state-resolved betas from legacy Gamma files.

17/07/24

"""

import pandas as pd
import numpy as np

def setTestMatE(gammaDF, rand=False, phase=True):
    """
    Set test matrix elements to Pandas DF.
    
    Matches (l,m,lam) terms from input gammas.
    
    Default case assigns 0.1*l+0.1
    
    If rand=True, assign random values
    
    If phase=True, assign random phase.
    
    """
    
    matE = []

    for l1 in gammaDF.index.levels[2]:
        for lam1 in gammaDF.index.levels[3]:
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