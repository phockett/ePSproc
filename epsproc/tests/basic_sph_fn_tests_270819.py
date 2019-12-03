# -*- coding: utf-8 -*-
# %%
"""
Script for testing spherical functions & products in ePSproc

27/08/19
"""

# %% Package imports
# %% Imports
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

# %% Function defs

# Define in root namespace
# runfile('E:/code/ePSproc/ePSproc/epsproc/ePSproc_IO.py', wdir='E:/code/ePSproc/ePSproc/epsproc')

# Define via import
import epsproc as ep

# %% Set files

# Set data dir
dataPath = os.path.join(os.getcwd(), 'data')

fList = ep.getFiles(fileBase = dataPath)

dataSet = ep.readMatEle(fileIn = fList[1])  # Load NO2 matrix elements for testing.

# %% Calculate spherical functions

# Spherical harmonics
YlmX = ep.sphCalc(6, res = 50)

# Wigner D
wDX = ep.wDcalc(eAngs = np.array([0,0,0]))

# %% Test direct multiplications

# Set matrix elements
daRed = ep.matEleSelector(dataSet[0], thres = 1e-2)

# Set Wigner D's
wDXre = wDX.sel({'lp':1, 'mu0':0})

testRe = wDXre.conjugate() * daRed

# Check terms by symmetry - for NO2 case this should select Cont = A2 case only
testRe.sel({'Cont':'A2'}).max()  # Should be non-zero

testRe.sel({'Cont':'B2'}).max() # Should be zero
testRe.sel({'Cont':'B1'}).max() # Should be zero

# %% Check MFPAD directly by multiplication...
#   MFPAD is matE * D * Ylm

testReY = testRe * YlmX
testPlot = testReY.sel({'Type':'L','it':1}).sum('mu') # Select on remaining dims
ep.sphSumPlotX(testPlot.squeeze().sum('Sym'), pType = 'p')

# %% Check (x,y,z) pol geometries

pRot = [0, 0, np.pi/2]
tRot = [0, np.pi/2, np.pi/2]
cRot = [0, 0, 0]
eAngs = np.array([pRot, tRot, cRot]) 

wDX = ep.wDcalc(eAngs = eAngs)
wDXre = wDX.sel({'lp':1, 'mu0':0})

testRe = wDXre.conjugate() * daRed * YlmX.conjugate()

testPlot = testRe.sel({'Type':'L','it':1}) # Select on remaining dims

for n in range(0,3):      
    # ep.sphSumPlotX(testPlot[n].squeeze().sum('Sym').sum('mu'), pType = 'a') # This runs for Euler as final dim, although not getting all correct results?
    ep.sphSumPlotX(testPlot[:,n].sum('Sym').sum('mu').sum('LM').squeeze(), pType = 'a')
