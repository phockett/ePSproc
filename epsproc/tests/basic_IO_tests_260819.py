# -*- coding: utf-8 -*-
"""
File for testing ePSproc basic IO

Created on Mon Aug 26 14:36:29 2019

@author: femto
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

file = fList[0]
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