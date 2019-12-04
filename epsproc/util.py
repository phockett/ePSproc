# -*- coding: utf-8 -*-
"""
ePSproc utility functions

Collection of small functions for sorting etc.

14/10/19    Added string replacement function (generic)
11/08/19    Added matEleSelector

"""

import numpy as np
import re
import scipy.constants

# Package fns.
from epsproc.basicPlotters import molPlot

#*************** Summary & display functions

# Print some jobInfo stuff & plot molecular structure
def jobSummary(jobInfo = None, molInfo = None, tolConv = 1e-2):
    """
    Print some jobInfo stuff & plot molecular structure. (Currently very basic.)

    Parameters
    ------------
    jobInfo : dict, default = None
        Dictionary of job data, as generated by :py:function:`epsproc.IO.headerFileParse()` from source ePS output file.

    molInfo : dict, default = None
        Dictionary of molecule data, as generated by :py:func:`epsproc.IO.molInfoParse()` from source ePS output file.

    tolConv : float, default = 1e-2
        Used to check for convergence in ExpOrb outputs, which defines single-center expansion of orbitals.

    Returns
    -------
    Job info : list

    """

    # Pull job summary info
    if jobInfo is not None:
        print('\n*** Job summary data')
        [print(line.strip('#')) for line in jobInfo['comments'][0:4]]
        print(f"\nElectronic structure input: {jobInfo['Convert'][0].split()[1].strip()}")
        print(f"Initial state occ:\t\t {jobInfo['OrbOccInit']}")
        print(f"Final state occ:\t\t {jobInfo['OrbOcc']}")
        print(f"IPot (input vertical IP, eV):\t\t {jobInfo['IPot']}")

        # Additional orb info
        print("\n*** Additional orbital info (SymProd)")
        iList = jobInfo['OrbOccInit'] - jobInfo['OrbOcc']
        print(f"Ionizing orb:\t\t\t {iList}")

        if molInfo is not None:
            # Get orbGrp for event
            iOrbGrp = np.flatnonzero(iList) + 1

            # Find entries in orbTableX
            orbSym = molInfo['orbTable'].where(molInfo['orbTable'].coords['OrbGrp'] == iOrbGrp, drop = True).coords['Sym']
            print(f"Ionizing orb sym:\t\t {np.unique(orbSym)}")

            orbIP = molInfo['orbTable'].where(molInfo['orbTable'].coords['OrbGrp'] == iOrbGrp, drop = True).sel(props = 'E')
            orbIPH = molInfo['orbTable'].where(molInfo['orbTable'].coords['OrbGrp'] == iOrbGrp, drop = True).sel(props = 'EH')
            print(f"Orb energy (eV):\t\t {np.unique(orbIP)}")
            print(f"Orb energy (H):\t\t\t {np.unique(orbIPH)}")

            orbIPnm = conv_ev_nm(np.unique(orbIP))
            print(f"Orb energy (cm^-1):\t\t {1/orbIPnm*1e7}")
            print(f"Threshold wavelength (nm):\t {np.abs(orbIPnm)[0]}")

            # Check ExpOrb outputs...
            ind = (molInfo['orbTable'][:,8].values < 1-tolConv) + (molInfo['orbTable'][:,8].values > 1+tolConv)
            if ind.any():
                print(f"\n*** Warning: some orbital convergences outside single-center expansion convergence tolerance ({tolConv}):")
                print(molInfo['orbTable'][ind, [0, 8]].values)



    # Display structure
    if molInfo is not None:
        print('\n*** Molecular structure\n')
        molPlot(molInfo)

    return jobInfo

# Print (LM) and symmetry sets with Pandas tables
def lmSymSummary(data):
    """Display summary info data tables.

    Works nicely in a notebook cell, with Pandas formatted table... but not from function?

    """

    print('\n*** Index summary\n')

    test = daPlot.Sym.to_pandas()
    print('Symmetry sets')
    print(test.unstack())

    test = daPlot.LM.to_pandas()
    print('\n(L,M) sets')
    print(test.unstack().T)




#*************** Selection functions

# Selector function for matrix elements in Xarray
def matEleSelector(da, thres = None, inds = None, dims = None, sq = False, drop=True):
    """
    Select & threshold raw matrix elements in an Xarray

    Parameters
    ----------
    da : Xarray
        Set of matrix elements to sub-select
    thres : float, optional, default None
        Threshold value for abs(matElement), keep only elements > thres.
        This is *element-wise*.
    inds : dict, optional, default None
        Dicitonary of additional selection criteria, in name:value format.
        These correspond to parameter dimensions in the Xarray structure.
        E.g. inds = {'Type':'L','Cont':'A2'}
    dims : str or list of strs, dimensions to look for max & threshold, default None
        Set for *dimension-wise* thresholding. If set, this is used *instead* of element-wise thresholding.
        List of dimensions, which will be checked vs. threshold for max value, according to abs(dim.max) > threshold
        This allows for consistent selection of continuous parameters over a dimension, by a threshold.
    sq : bool, optional, default False
        Squeeze output singleton dimensions.
    drop : bool, optional, default True
        Passed to da.where() for thresholding, drop coord labels for values below threshold.

    Returns
    -------
    daOut
        Xarray structure of selected matrix elements.
        Note that Nans are dropped if possible.

    Example
    -------
    >>> daOut = matEleSelector(da, inds = {'Type':'L','Cont':'A2'})

    """

    # Iterate over other selection criteria
    # This may return view or copy - TBC - but seems to work as expected.
    # http://xarray.pydata.org/en/v0.12.3/indexing.html#copies-vs-views
    if inds is not None:
        da = da.sel(inds)    # Fors inds as dict, e.g. {'Type':'L','it':1,'Cont':'A2'}
                                # May want to send as list, or automate vs. dim names?
                                # NOTE - in current dev code this is used to reindex, so .squeeze() casuses issues!


    # Reduce dims by thesholding on abs values
    # Do this after selection to ensure Nans removed.
    if thres is not None and dims is None:
        daOut = da.where(np.abs(da) > thres, drop = drop)
    else:
        daOut = da

    # If dims is set, check over dims for consistency.
    # WILL this just produce same results as thres then squeeze...?
    if dims is not None:
        daOut = daOut.where(np.abs(da).max(dim = dims) > thres, drop = drop)

    if sq:
        daOut = daOut.squeeze()  # Squeeze dims.

    return daOut


# Select over vals from data structure (list)
# Currently only used in IO.matEleGroupDim
def dataGroupSel(data, dInd):
    a = data[0]
    dataSub = []

    uVals = np.unique(a[dInd,:])

    for val in uVals:
        # Get matching terms and subset data
        # iSel = np.nonzero(a[dInd,:]==val)
        iSel = (a[dInd,:]==val)
        dataSub.append([data[0][:,iSel], data[1][iSel]])

    return dataSub


# Xarray groupby + compare values
# STARTED... but not finished.  For basic diff along a dimension, just use da.diff(dim), see http://xarray.pydata.org/en/stable/generated/xarray.DataArray.diff.html#xarray.DataArray.diff
# def groupCmp(data, dim):
#     """
#     Basic routine to compare sets of values by dimension, using Xarray groupby functionality.
#
#     Parameters
#     ----------
#     data : Xarray
#         Data for comparison
#
#     dim : str
#         Dimension label for grouping
#
#     Returns
#     -------
#
#     """
#
#     dGroup = data.groupby(dim)
#
#     # Check differences between groups
#     for gTest in dGroup:



#************************** LIST functions

# Return list of standard dataArray dims for matrix elements
def matEdimList(sType = 'stacked'):
    """
    Return standard list of dimensions for matrix elements.

    Parameters
    ----------
    sType : string, optional, default = 'stacked'
        Selected 'stacked' or 'unstacked' dimensions.
        Set 'sDict' to return a dictionary of unstacked <> stacked dims mappings for use with xr.stack({dim mapping}).

    Returns
    -------
    list : set of dimension labels.

    """
    if sType is 'stacked':
        # stackedDims
        return ['LM', 'Eke', 'Sym', 'mu', 'it', 'Type']

    elif sType is 'sDict':
        return {'LM':['l','m'],'Sym':['Cont', 'Targ', 'Total']}

    else:
        # unstackedDims
        return ['l','m', 'Eke', 'Cont', 'Targ', 'Total', 'mu', 'it', 'Type']

# Return list of standard dataArray dims for BLM values
def BLMdimList(sType = 'stacked'):
    """
    Return standard list of dimensions for calculated BLM.

    Parameters
    ----------
    sType : string, optional, default = 'stacked'
        Selected 'stacked' or 'unstacked' dimensions.
        Set 'sDict' to return a dictionary of unstacked <> stacked dims mappings for use with xr.stack({dim mapping}).

    Returns
    -------
    list : set of dimension labels.

    """
    if sType is 'stacked':
        # stackedDims
        return ['Euler', 'Eke', 'BLM']

    elif sType is 'sDict':
        return {'BLM':['l','m'],'Euler':['P','T','C']}

    else:
        # unstackedDims
        return ['Eke', 'l', 'm', 'P', 'T', 'C']

# Return list of standard dataArray dims for Euer angles
def eulerDimList(sType = 'stacked'):
    """
    Return standard list of dimensions for frame definitions, from :py:func:`epsproc.sphCalc.setPolGeoms()`.

    Parameters
    ----------
    sType : string, optional, default = 'stacked'
        Selected 'stacked' or 'unstacked' dimensions.
        Set 'sDict' to return a dictionary of unstacked <> stacked dims mappings for use with xr.stack({dim mapping}).

    Returns
    -------
    list : set of dimension labels.

    """
    if sType is 'stacked':
        # stackedDims
        return ['Euler']

    elif sType is 'sDict':
        return {'Euler':eulerDimList(sType = 'unstacked')}

    else:
        # unstackedDims
        return ['Label','P','T','C']

# Return list of standard dataArray dims for Euer angles
def ADMdimList(sType = 'stacked'):
    """
    Return standard list of dimensions for frame definitions, from :py:func:`epsproc.sphCalc.setADMs()`.

    Parameters
    ----------
    sType : string, optional, default = 'stacked'
        Selected 'stacked' or 'unstacked' dimensions.
        Set 'sDict' to return a dictionary of unstacked <> stacked dims mappings for use with xr.stack({dim mapping}).

    Returns
    -------
    list : set of dimension labels.

    """
    if sType is 'stacked':
        # stackedDims
        return ['ADM','t']

    elif sType is 'sDict':
        return {'ADM':['K','Q','S'],'t':'t'}

    else:
        # unstackedDims
        return ['K','Q','S','t']




# Return a dict of allowed dataTypes, corresponding to ePS commands/file segments, or epsproc processed data.
def dataTypesList():
    """
    Return a dict of allowed dataTypes, corresponding to epsproc processed data.

    Each dataType lists 'source', 'desc' and 'recordType' fields.

    - 'source' fields correspond to ePS functions which get or generate the data.
    - 'desc' brief description of the dataType.
    - 'recordType' gives the required segment in ePS files (and associated parser). If the segment is not present in the source file, then the dataType will not be available.

    TODO: best choice of data structure here?  Currently nested dictionary.

    """

    dataDict = {'BLM' :
                    {'source':'epsproc.MFBLM',
                    'desc':'Calcualted MF beta parameters from epsproc.MFBLM(), based on dipole matrix elements from ePS.',
                    'recordType':'DumpIdy',
                    'dims': BLMdimList(sType = 'sDict')
                    },
                'matE' :
                    {'source':'epsproc.IO.readMatEle',
                    'desc':'Raw photoionization matrix elements from ePS, DumpIdy command and file segments.',
                    'recordType':'DumpIdy',
                    'dims': matEdimList(sType = 'sDict')
                    },
                'EDCS' :
                    {'source':'epsproc.IO.readMatEle',
                    'desc':'Electron scattering DCS results from ePS, EDCS command and file segments.',
                    'recordType':'EDCS'
                    },
                'XSect' :
                    {'source':'CrossSection',
                    'desc':'Photoionziation cross section results from ePS, GetCro command and CrossSection file segments.',
                    'recordType':'CrossSection'
                    },
                'Euler' :
                    {'source':'epsproc.sphCalc.setPolGeoms()',
                     'desc':'Frame rotation definitions in terms of Euler angles and corresponding Quaternions.',
                     'dims':eulerDimList(sType = 'sDict')
                     },
                 'ADM' :
                     {'source':'epsproc.sphCalc.setADMs()',
                      'desc':'Defines ADMs, A(K,Q,S), for aligned frame distributions.',
                      'dims':ADMdimList(sType = 'sDict')
                      }
                }

    return dataDict


#***************** Convenience functions...

# Multistring replace
# See https://gist.github.com/bgusach/a967e0587d6e01e889fd1d776c5f3729
# https://stackoverflow.com/questions/6116978/how-to-replace-multiple-substrings-of-a-string
def stringRepMap(string, replacements, ignore_case=False):
    """
    Given a string and a replacement map, it returns the replaced string.
    :param str string: string to execute replacements on
    :param dict replacements: replacement dictionary {value to find: value to replace}
    :param bool ignore_case: whether the match should be case insensitive
    :rtype: str

    CODE from:
    https://gist.github.com/bgusach/a967e0587d6e01e889fd1d776c5f3729
    https://stackoverflow.com/questions/6116978/how-to-replace-multiple-substrings-of-a-string
    ... more or less verbatim.

    Thanks to *bgusach* for the Gist.

    """
    # If case insensitive, we need to normalize the old string so that later a replacement
    # can be found. For instance with {"HEY": "lol"} we should match and find a replacement for "hey",
    # "HEY", "hEy", etc.
    if ignore_case:
        def normalize_old(s):
            return s.lower()

        re_mode = re.IGNORECASE

    else:
        def normalize_old(s):
            return s

        re_mode = 0

    replacements = {normalize_old(key): val for key, val in replacements.items()}

    # Place longer ones first to keep shorter substrings from matching where the longer ones should take place
    # For instance given the replacements {'ab': 'AB', 'abc': 'ABC'} against the string 'hey abc', it should produce
    # 'hey ABC' and not 'hey ABc'
    rep_sorted = sorted(replacements, key=len, reverse=True)
    rep_escaped = map(re.escape, rep_sorted)

    # Create a big OR regex that matches any of the substrings to replace
    pattern = re.compile("|".join(rep_escaped), re_mode)

    # For each match, look up the new string in the replacements, being the key the normalized old string
    return pattern.sub(lambda match: replacements[normalize_old(match.group(0))], string)


#********************** Calculations

# Convert eV <> Hartrees (atomic units)
def conv_ev_atm(data, to = 'ev'):
    """
    Convert eV <> Hartree (atomic units)

    Parameters
    ----------
    data : int, float, np.array
        Values to convert.

    to : str, default = 'ev'
        - 'ev' to convert H > eV
        - 'H' to convert eV > H

    Returns
    -------
    data converted in converted units.

    """

    # Define H in eV, value from https://en.wikipedia.org/wiki/Hartree_atomic_units#Units
    # H = 27.211386245
    H = scipy.constants.physical_constants['Hartree energy in eV'][0]  # Use scipy value

    if to is 'ev':
        dataOut = data*H
    else:
        dataOut = data/H

    return dataOut

# Convert energy in eV to wavelength in nm
def conv_ev_nm(data): #, to = 'nm'):
    """Convert E(eV) <> nu(nm)."""

    # Define constants from scipy.constants
    h = scipy.constants.h
    c = scipy.constants.c
    evJ = scipy.constants.physical_constants['electron volt-joule relationship'][0]

    # Define output units - wavelength in m
    waveConv = 1e-9
    dataOut = (h * c)/(data * evJ)/waveConv

    # if to is 'nm':
    #     dataOut = (h * c)/(data * evJ)/waveConv
    #
    # else:
    #     dataOut = (data/waveConv)*evJ/(h * c)

    return dataOut
