"""
ePSproc Gamma functions: quick plotters

NOTE: these assume Pandas DataFrames as output by gammaCalc functions.

For more general plotters, see epsproc.basicPlotters and class-based plotters.

TODO: integrate PD and XR plotters...?

09/08/24

"""


def betaCalcPlot(betaDF, colName='t', thres=None):
    """
    Quick hvplot() for Pandas DataFrame outputs from gammaCalc.betaCalc().
    
    Plots (L,M) values vs. index (default 't').
    
    """
    
    if thres is not None:
        # Threshold and tidy
        betaThres = betaDF.where(np.abs(betaDF.max(axis=1))>1e-6).dropna()
    else:
        betaThres = betaDF

    # Reformat for plotting
    betaThres.columns.name = 't'
    betaPlot = betaThres.apply(np.real)   # TODO: use normal plotSelector here.
    
    plotT = betaPlot.reset_index()
    plotT.set_index(plotT['L'].astype(str)+','+plotT['M'].astype(str),inplace=True)
    plotT = plotT.T
    
    # Plot and return - note skip first two rows which are here (L,M) values.
    return plotT[2:].hvplot.line()