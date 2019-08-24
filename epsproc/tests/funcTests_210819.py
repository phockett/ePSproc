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

#%% Subselections using matEleSelector
#   THIS IS UGLY, but seems to work - get da of correct dims out.


indList = ['ip','it','mu']

da = dataSet[0].copy()

daRedList = []
for x in np.unique(da[indList[0]]):
    daRedList0 = []
    for y in np.unique(da[indList[1]]):
        daRedList1 = []
        for z in np.unique(da[indList[2]]):
            red = ep.UT.matEleSelector(da, inds = {indList[0]:x, indList[1]:y, indList[2]:z})
            
            red = red.expand_dims({indList[0]:[x], indList[1]:[y], indList[2]:[z]})
            
            daRedList1.append(red)
            
        daOut1 = xr.combine_nested(daRedList1, concat_dim = indList[2])
        daRedList0.append(daOut1)    
        
    daOut2 = xr.combine_nested(daRedList0, concat_dim = indList[1])
    daRedList.append(daOut2)

daOut =  xr.combine_nested(daRedList, concat_dim = indList[0])
    
#%% Test dataset

daOut2 = daOut.copy()


test = xr.Dataset(daOut, daOut2)
        