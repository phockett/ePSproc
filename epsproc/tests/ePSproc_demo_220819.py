# -*- coding: utf-8 -*-
"""
ePSproc basic demo/test script

TODO: Move to Jupyter notebook for distro

22/08/19

"""
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

#%% ePSproc

import epsproc as ep

# Set data dir
dataPath = os.path.join(os.getcwd(), 'data')

#%% Read data

dataSet = ep.readMatEle(fileBase = dataPath)  # Scan data dir

#%% Calculate MFPADs

TaX, TlmX = ep.mfpad(dataSet[0])

#%% Plot

# Plot for each set of Euler angles, assuing they are top-level index
# Should be a neater way to do this...?
for n in range(0,3):
    ep.sphSumPlotX(TaX[n].sum(dim = 'ES'))
