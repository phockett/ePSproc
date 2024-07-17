"""
ePSproc Gamma functions: IO

- Compute state-resolved betas from legacy Gamma files.

17/07/24

"""

import pandas as pd

def setTestMatE(gammaDF):
    """
    Set test matrix elements to Pandas DF.
    
    Matches (l,m,lam) terms from input gammas.
    
    """
    
    matE = []

    for l1 in gammaDF.index.levels[2]:
        for lam1 in gammaDF.index.levels[3]:
            matE.append([l1,lam1, 0.1+l1*0.1])

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
    


def assignMatE(gammaDF, matE=None):
    """
    Assign matrix elements as columns in Pandas DataFrame of gamma values.
    
    """
    
    if matE is None:
        matE1 = setTestMatE(gammaDF)
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

    return df2



def betaCalc(gammaDF, matE=None):
    """
    Compute betas from (legacy) gamma terms and matrix elements.
    
    All values must be as Pandas DataFrames
    
    """
    
    # Assign all terms to master DF
    dfCalc = assignMatE(gammaDF,matE)
    
    # Compute product terms
    # Multiply cols
    dfCalc['product'] = dfCalc['betaTerm']*dfCalc['gamma']*dfCalc['matE1']*dfCalc['matE2']
    
    # Sum terms
    dfSum = dfCalc.groupby(by=['L','M']).sum()
    
    return dfSum, dfCalc