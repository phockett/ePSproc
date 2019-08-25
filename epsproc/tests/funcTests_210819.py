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
#   THIS IS UGLY, but seems to work - get da of correct dims out (with multilevel coords in).
#   Could also try dataset to array for split and recombine...?
#   http://xarray.pydata.org/en/v0.12.3/reshaping.html


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

#%% Tidier version (maybe) with reduced nesting
#   Linear tree ip > mu >it

# def matEleSplit(da, indList)
daIn = dataSet[0].copy()
daRedList = []


# Split on 'ip' - will always be (1,2), and split matEle into two
ipLabel = ['L','V']
for n, val in enumerate(range(1,3)):
    tmp = ep.UT.matEleSelector(daIn, inds = {'ip':val})
    tmp = tmp.expand_dims({'Type':[ipLabel[n]]})
    daRedList.append(tmp)

# Restack
daRed = xr.combine_nested(daRedList, concat_dim = 'Type')

# Split on mu - values from set {-1,0,1} depending on symmetry
daRedList = []
uVals = np.unique(daRed.mu)
for n, val in enumerate(uVals):
    tmp = ep.UT.matEleSelector(daRed, inds = {'mu':val})
    tmp = tmp.expand_dims({'mu':[val]})
    daRedList.append(tmp)

# Restack
daRed = xr.combine_nested(daRedList, concat_dim = 'mu')    

# Split on it
daRedList = []
uVals = np.unique(daRed.it)
for n, val in enumerate(uVals):
    tmp = ep.UT.matEleSelector(daRed, inds = {'it':val})
    tmp = tmp.expand_dims({'it':[val]})
    daRedList.append(tmp)

# Restack
daRed = xr.combine_nested(daRedList, concat_dim = 'it')    



#%% Using datasets
#   Should work... but notation gets ugly(er)
#   Actually - used arrays for all coords which will be used in a sum.
#   The datasets to hold remaining dims (ip)...?
#   WAIT -- in this construction will still require knowldege of dims.

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
            
        # daOut1 = xr.Dataset({'mu0':daRedList1[0], 'mu1':daRedList1[1]})
        daOut1 = xr.combine_nested(daRedList1, concat_dim = indList[2])
        daRedList0.append(daOut1)    
        
    daOut2 = xr.combine_nested(daRedList0, concat_dim = indList[1])
    daRedList.append(daOut2)

daOut =  xr.combine_nested(daRedList, concat_dim = indList[0])
# daOut = xr.Dataset({})

#%% Recursive form...

# Define multi-level selector to recursively parse inds and concat dims
# ALMOST WORKS... issue is what to pass recursively, once there are multiple arrays created.
def matEleSelectorRec(da, inds):
    
    indRed = inds.copy()
    level = indRed.pop()
    daRedList = []
    
    if len(indRed)>0:
        daRedList.append(matEleSelectorRec(da, indRed))
        
    else:
        for daRedIn in daRedList:
            for val in np.unique(daRedIn[level]):
                daRed = ep.UT.matEleSelector(daRedIn, inds = {level:val})
                daRed = daRed.expand_dims({level:[val]})
                
                daRedList.append(daRed)

    return daRedList

indList = ['ip','it','mu']

da = dataSet[0].copy()

daRedList = matEleSelectorRec(da,indList)


#%% Test E reorg
tmp = daOut.sel({'E':0.81})
tmp2 = tmp.copy()
tmp = tmp.expand_dims({'Eke':[1]})
tmp2 = tmp2.expand_dims({'Eke':[2]})

daE = xr.combine_nested([tmp,tmp2], concat_dim = 'Eke')


#%% Test dataset

ds1 = xr.Dataset({'matE':daE})

ds2 = xr.Dataset({'E1':tmp,'E2':tmp2})

# Seems to autostack - lists all Eke here.
# May be able to use for sorting above...?
ds2['E2']