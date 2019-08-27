# -*- coding: utf-8 -*-
"""
Script for testing ePSproc basic IO

26/08/19

"""

#%% Package imports
#%% Imports
import os
import re
import numpy as np
import pandas as pd
from io import StringIO
import xarray as xr
from scipy.special import sph_harm

# For plotting functions
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

from matplotlib import cm, colors

#%% Function defs

# Define in root namespace
# runfile('E:/code/ePSproc/ePSproc/epsproc/ePSproc_IO.py', wdir='E:/code/ePSproc/ePSproc/epsproc')

# Define via import
import epsproc as ep

#%% Set files

# Set data dir
dataPath = os.path.join(os.getcwd(), 'data')

fList = ep.getFiles(fileBase = dataPath)

#%% Read file

file = fList[1]  # HERE SET FOR [0] for N2 multi-E data, [1] for NO2 single E test set.
lines, dumpSegs = ep.dumpIdyFileParse(file)

#%% Check parameters

ekeList = ep.scatEngFileParse(file)

symSegs = ep.symFileParse(file)

#%% Parse matrix elements - np array format

matEleList = []
attribsList = []

for dumpSeg in dumpSegs:
    matEle, attribs = ep.dumpIdySegParse(dumpSeg)
    matEleList.append(matEle)
    attribsList.append(attribs)


#%% Parse matrix elements - Xarray format (already includes loop over segments)

data, blankSegs = ep.dumpIdySegsParseX(dumpSegs, ekeList, symSegs)

#%% Full file IO

dataSet = ep.readMatEle(fileIn = file)

#%% Test sph multiplication over all dims

res = 50
# daRed = ep.matEleSelector(dataSet[0], thres = 1e-2, inds = {'Type':'L','it':1})
daRed = ep.matEleSelector(dataSet[0], thres = 1e-2, inds = {})

# Generate spherical harmonics
Lmax = daRed.l.max()
YlmX = ep.sphCalc(Lmax, res = res)

# Multiply
test = daRed * YlmX

# Check some terms by selecting a single set of matrix elements...
inds = {'mu':0, 'Eke':0.1, 'Cont':'SU'}
daRed.sel(inds)

# Values should match file:
#   m   l mu ip it Value
#   0   1  0  1  1  0.27363206E+01 -0.92696864E-01
#   0   3  0  1  1 -0.17158793E+00 -0.79529163E+00

test.sel(inds)
ep.sphSumPlotX(test.sel(inds).squeeze(dim = 'Sym'))

#%% Test Xarray image plot with faceting
# http://xarray.pydata.org/en/stable/plotting.html#dimensional

# Facet by E
inds = {'mu':0, 'Type':'L', 'Cont':'SU'}
testPlot = np.abs(test.sel(inds).squeeze(dim = 'Sym').squeeze(dim = 'it').sum(dim='LM'))
testPlot = testPlot.isel(Eke=slice(0, 50, 5))   # Slice data
testFacet = testPlot.plot(x='Theta', y='Phi', col='Eke', col_wrap=5)

# Facet by (E, type)
inds = {'mu':0, 'Cont':'SU'}
testPlot = np.abs(test.sel(inds).squeeze(dim = 'Sym').squeeze(dim = 'it').sum(dim='LM'))
testPlot = testPlot.isel(Eke=slice(0, 50, 10))   # Slice data
testFacet = testPlot.plot(x='Theta', y='Phi', col='Eke', row='Type', col_wrap=5)