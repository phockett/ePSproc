# -*- coding: utf-8 -*-
"""
Script for testing ePSproc packaged functions.

21/08/19

"""

#%% Import

import epsproc as ep

#%% Read (from /data)

# dataSet = ep.IO.readMatEle()
dataSet = ep.readMatEle()


#%% MFPADs

#ep.MF.mfpad(dataSet[0])
ep.mfpad(dataSet[0])