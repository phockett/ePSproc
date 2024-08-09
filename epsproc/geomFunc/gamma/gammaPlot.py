"""
ePSproc Gamma functions: quick plotters

NOTE: these assume Pandas DataFrames as output by gammaCalc functions.

For more general plotters, see epsproc.basicPlotters and class-based plotters.

TODO: integrate PD and XR plotters...?

09/08/24

"""

import numpy as np
from epsproc.sphPlot import plotTypeSelector

def betaCalcPlot(betaDF, colName='t', thres=None, pType='a'):
    """
    Quick hvplot() for Pandas DataFrame outputs from gammaCalc.betaCalc().
    
    Plots (L,M) values vs. index (default 't').
    
    TODO: integrate PD and existing XR plotters...?
    
    
    Parameters
    ----------
    colName : str, optional, default='t'
        Set name for columns.
        Default case has cols per t.
        
    thres : float, optional, default = None
        If set, threshold abs values.
        
    pType : str, optional, default = 'a'
        Data type to plot, per :py:func:`epsproc.sphPlot.plotTypeSelector()`.
        Default is abs values.
    
    
    Return
    ------
    hvplot object
    
    """
    
    if thres is not None:
        # Threshold and tidy
        betaThres = betaDF.where(np.abs(betaDF.max(axis=1))>thres).dropna()
    else:
        betaThres = betaDF

    # Reformat for plotting
    betaThres.columns.name = colName
    
    # betaPlot = betaThres.apply(np.real)   # TODO: use normal plotSelector here.
    # Set plot type...
    pTypeDict = plotTypeSelector(returnDict=True)
    
    if pType in pTypeDict.keys():
        betaThres = betaThres.apply(pTypeDict[pType]['Exec'])
    else:
        print(f"*** pType = {pType} not implement, defaulting to 'a' (abs values). See epsproc.plotTypeSelector for options.")
        betaThres = betaThres.apply(pTypeDict['a']['Exec'])
    
    
    plotT = betaThres.reset_index()
    plotT.set_index(plotT['L'].astype(str)+','+plotT['M'].astype(str),inplace=True)
    plotT = plotT.T
    
    # Plot and return - note skip first two rows which are here (L,M) values.
    return plotT[2:].hvplot.line()