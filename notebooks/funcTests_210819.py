# -*- coding: utf-8 -*-
"""
Script for testing ePSproc packaged functions.

21/08/19

"""

#%% Import
import matplotlib.pyplot as plt
import epsproc as ep

#%% Read (from /data)

# dataSet = ep.IO.readMatEle()
dataSet = ep.readMatEle()


#%% MFPADs

#ep.MF.mfpad(dataSet[0])

TX, TlmX = ep.mfpad(dataSet[1])

# Plot for each pol geom (symmetry)
for n in range(0,3):      
    # ep.sphSumPlotX(testPlot[n].squeeze().sum('Sym').sum('mu'), pType = 'a') # This runs for Euler as final dim, although not getting all correct results?
    ep.sphSumPlotX(TX[n].sum('Sym').squeeze(), pType = 'a')

# Plot abs(TX) images using Xarray functionality
TX.squeeze().pipe(np.abs).plot(x='Theta',y='Phi', col='Euler', row='Sym')

#%% Plot matrix elements using Xarray functionality

# daRed = ep.matElementSelector(dataSet[1], inds = {})
# daPlot = np.abs(dataSet[0].sum('mu').sum('Sym').sel({'Type':'L'}).squeeze())
daPlot = dataSet[0].sum('mu').sum('Sym').sel({'Type':'L'}).squeeze()

fig, axes = plt.subplots(nrows=2)

daPlot.pipe(np.abs).plot.line(x='Eke', ax=axes[0])

# Doesn't work for angle - outputs np.array
# daPlot.pipe(np.angle).plot.line(x='Eke', ax=axes[1])

#%% As colourmap - doens't work for category data?

daPlot.plot.pcolormesh('LM', 'Eke')

#%% Plot with faceting on type

daPlot = dataSet[0].sum('mu').sum('Sym').squeeze()

daPlot.pipe(np.abs).plot.line(x='Eke', row='Type')

#%% Plot with faceting on symmetry

daPlot = dataSet[0].sum('mu').squeeze()

daPlot.pipe(np.abs).plot.line(x='Eke', col='Sym', row='Type')

#%% Plot MFPAD surfaces vs E

TX, TlmX = ep.mfpad(dataSet[0])

TXplot = TX.sum('Sym').squeeze().isel(Eke=slice(0,50,10))
TXplot.pipe(np.abs).plot(x='Theta',y='Phi', row='Euler', col='Eke')

