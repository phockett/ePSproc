# -*- coding: utf-8 -*-
"""
ePSproc convenience functions

Collection of small functions for sorting etc.


"""

# import numpy as np
import re
import itertools
import os
from datetime import datetime

# import itertools   # For itertools.chain.from_iterable
from collections.abc import Iterable   # For iterable checking

# import scipy.constants
#
# # Package fns.
# from epsproc.basicPlotters import molPlot

# Used for data type checking later, but may not be required if this can be done per data object?
import xarray as xr
import pandas as pd

from epsproc.util.listFuncs import getRefDims


try:
    from natsort import natsorted  # For natural sorting
    natsortFlag = True

except ImportError as e:
    if e.msg != "No module named 'natsort'":
        raise
    print('* natsort not found, some sorting functions not available. ')
    natsortFlag = False

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

# Sort a 2D numpy array.
def arraySort2D(a, col):
    """
    Sort np.array `a` by specified column `col`.
    From https://thispointer.com/sorting-2d-numpy-array-by-column-or-row-in-python/
    """
    return a[a[:,col].argsort()]


# Set up lambdas for itertools groupby use in fileListSort (below)
# Quick hack to get this working for different file-naming conventions.
# TODO: make nicer & generalise.
# TODO: consider cases with/without prefix str for single and multiple dirs - that's the main difference with prefixStr...?
# 28/04/21 - currently broken for wavefn files, must have changed this for other purposes AFTER https://epsproc.readthedocs.io/en/dev/demos/ePSproc_wfPlot_tests_150720-110820-CH3I-tidy_Stimpy.html ?
#           Quick fix by also matching by file type for orb data files (.dat)
#           Should really use regex here!
def sortGroupFn(fListSorted, prefixStr):

    # (1) Original case, works for wavefunction files with naming convention
    #  <jobSym>_<Eke>_Orb.dat, e.g. CH3ISA1CA1_1.0eV_Orb.dat
    #  In this case, split and group on first part of file name
    partName = fListSorted[0].replace(prefixStr,'')
    if (len(partName.split('_')) < 2) or (partName.endswith('.dat')):
        return lambda x:x.replace(prefixStr,'').split('_')[0]

    # (2) Case for multi-E ePS job output files.
    #  <job>.<orb>_<Sym>_<Eke>.out, e.g. N2O_wf.orb1_S_E1.0_6.0_97.0eV.inp.out
    # In this case, just group be prefix, which should be OK if only a single dir is set.
    # Should likely also check for file extension or other here?
    else:
        # return lambda x:prefixStr  # Use prefix str only
        return lambda x:x.split('_E')[0]  # Check from full name, no additional prefixStr required.



# Sort & group filenames
def fileListSort(fList, groupByPrefix=True, prefixStr = None, verbose=1):
    """
    Sort a list of file names, and group by prefix.

    Note: this currently assumes a file name schema whereby split('_')[0] picks the grouping string.

    Note: os.path.commonprefix() is used for determining prefix, this may fail in some cases (e.g. for cases where a single file is passed, or files from different dirs).
    Pass prefix manaully in these cases.

    Returns
    -------
    fListSorted, groupedList, prefixStr


    """

    if natsortFlag:
        fListSorted = natsorted(fList)
    else:
        fListSorted = sorted(fList)

    # prefixStr = ''
    if groupByPrefix:
        if prefixStr is None:
            prefixStr = os.path.commonprefix(fListSorted)  # Find common prefix if not passed.

        # Solution with itertools groupby
        # Adapted from https://stackoverflow.com/a/13368753
        # groupedList = [list(v) for k,v in itertools.groupby(fListSorted,key=lambda x:x.replace(prefixStr,'').split('_')[0])]
        groupedList = [list(v) for k,v in itertools.groupby(fListSorted,key=sortGroupFn(fListSorted, prefixStr))]


    if verbose:
        print(f"\n*** FileListSort \n  Prefix: {prefixStr} \n  {len(groupedList)} groups.")

        if verbose > 1:
            print("\n  Grouped list:")
            print(*groupedList, sep = '\n')

    if len(fList) > 1:
        return fListSorted, groupedList, prefixStr
    else:
        return fList, fList, None


# Return a time-string for setting unique file names
# May already have this elsewhere...?
def timeStamp():
    """Get local time and return formatted string "%d-%m-%y_%H-%M-%S" for time-stamping filesnames."""

    dt = datetime.now()

    return dt.strftime("%d-%m-%y_%H-%M-%S")



def checkDims(data, refDims = []):
    """
    Check dimensions for a data array (Xarray) vs. a reference list (or dict).

    Parameters
    ----------
    data : Xarray
        Data array to check.

    refDims : str, list, dict, optional, default = []
        Dims to check vs. input data array.
        If dict is passed only keys (==stacked dims) are tested.
        Update 06/06/22: now also checks unstacked case & returns safe ref mapping.

    Returns
    -------

    dictionary
        Containing:

        - stacked and unstacked dims
        - stacked dim mappings
        - intersection and differences vs. refDims
        - safeStack for use with restack() function even if dims are missing.

    Examples
    --------

    >>> # Return dim lists
    >>> ep.util.misc.checkDims(dataTest)

    >>> # Return dim lists
    >>> ep.util.misc.checkDims(dataTest)


    TODO: check and order dims by size? Otherwise set return is alphebetical

    06/06/22 Added better support for stacked refDims & remapping with refDims.
             Now also tests refDims items (unstacked dims), as well as keys (stacked dims).
             Added outputs 'extraUS', 'invalidUS', 'remap'
             May now be some repetition here, but didn't touch existing items to ensure back compatibility!

    28/09/21 Added basic support for Pandas DataFrames, still needs some work.
             See https://stackoverflow.com/questions/21081042/detect-whether-a-dataframe-has-a-multiindex for some thoughts.

    26/08/21 Added additional tests for stacked dims vs. ref.
             Bit messy, but left as-is to avoid breaking old code.
             In future should amalgamate stacked & unstacked tests & tidy output.

    23/08/21 Added stacked dim mapping output.
    11/05/21 Added handling for stacked dims.

    """

    # refDimsUS = []

    # Force to list to avoid breaking *unpacking later
    if not isinstance(refDims, (list, dict)):
        refDims = [refDims]
        # refDimsUS = []

    # For stacked case, also get unstacked dims.
    if isinstance(refDims, dict):
        # refDimsUS = list(itertools.chain.from_iterable(refDims.values()))  # This can fail for mixed types with non-iterables

        # More robust flattening for all types - probably there is a cleaner way to do this however?
        refDimsUS = []
        for v in refDims.values():
            if isinstance(v, Iterable):
                refDimsUS.extend(v)
            else:
                refDimsUS.append(v)

    else:
        refDimsUS = refDims  # Case for unstacked refDims


    # Get dim names from Xarray
    if isinstance(data, xr.DataArray):
        dims = data.dims # Set dim list - this excludes stacked dims
        dimsUS = data.unstack().dims  # Set unstaked (full) dim list

        stackedDims = list(set(dims) - set(dimsUS))

        stackedDimsMap = {k:data.indexes[k].names for k in stackedDims}  # Get stacked dim mapping from indexes (Xarray/Pandas)

    # Get BOTH columns and index names & items from Pandas DataFrame
    if isinstance(data, pd.DataFrame):
        # Basic version
        # if isinstance(data.index, pd.MultiIndex):
        #     print("*** Warning: checkDims doesn't yet support Pandas MultiIndex.")
        #     # TODO: add full MultiIndex unstack/pull names here. See https://stackoverflow.com/questions/21081042/detect-whether-a-dataframe-has-a-multiindex
        #     # .index.names should work?  But only for MultiIndex, or returns none? And for cols too?
        #     # UPDATE: now added more general routine below, but doesn't differentiate between stacked & unstacked dims.
        #     dims = []
        #     dimUS = []
        #
        # else:
        #     dims = data.columns.append(data.index.names) # Set dim list, assume this is both columns & rows
        #     dimUS = []

        # Get names & items from Pandas DataFrame
        # Assume that dim labels must be strings, skip if not, or if None
        # This currently gets both stacked & unstacked dims and puts these in a single list.
        dims = []
        dimsUS = []
        for item in ['columns','index']:
            base = getattr(data, item)  #.to_list()

            if isinstance(base, pd.MultiIndex):
                # print(f"*** Warning: found MultiIndex for DataFrame data.{item} - checkDims doesn't yet support Pandas MultiIndex.")
                print(f"*** Warning: found MultiIndex for DataFrame data.{item} - checkDims may have issues with Pandas MultiIndex, but will try anyway.")

            # Get data.item and data.item.names
            for nameList in [base.to_list(), list(getattr(base, 'names'))]:
                dimsUS.extend((nameList if ((nameList[0] is not None) and (isinstance(nameList[0], str))) else []))

        dims = dimsUS

        stackedDims = [] #list(set(dims) - set(dimsUS))

        stackedDimsMap = {} # {k:data.index[k].names for k in stackedDims}  # Get stacked dim mapping from indexes (Xarray/Pandas)


    # Check ref vs. full dim list (unstacked dims)
    sharedDims = list(set(dimsUS)&{*refDims})  # Intersection
    extraDims = list(set(dimsUS) - {*refDims})  # Difference
    invalidDims = list({*refDims} - set(dimsUS))  # Missing

    # Check also unstacked dims - this is useful in remapping cases
    extraDimsUS = list(set(dimsUS) - {*refDimsUS})  # Difference Unstacked
    invalidDimsUS = list({*refDimsUS} - set(dimsUS))  # Missing

    # 26/08/21 - added additional tests for stacked dims vs. ref, but may want to amalgamate with above?
    # TODO: tidy output too - separate dicts for stacked and unstacked dims?
    # Check ref vs. full dim list (stacked dims), note also unstacked dims subtracted to avoid duplicated dims.
    sharedDimsStacked = list(set(dims)&{*refDims} - set(dimsUS))  # Intersection
    extraDimsStacked = list(set(dims) - {*refDims} - set(dimsUS))  # Difference
    invalidDimsStacked = list({*refDims} - set(dims) - set(dimsUS))  # Missing

    # Test also missing dims overall
    missingDims = list({*refDims} - set(dimsUS) - set(dims))  # Missing

    # Set remapping dict for safe dims only
    # safeRemap = list({*refDims} - {*invalidDimsStacked})
    safeStack = {}

    # This works
    # for k,v in refDims.items():
    #     if k in invalidDimsStacked:
    #         safeRemap[k] = [item for item in v if item not in invalidDimsUS]

    # Better (?) - only check missing stacked dims
    for k in invalidDimsStacked:
        # if k in refDims.keys():
         stackList = [item for item in refDims[k] if item not in invalidDimsUS]

         if stackList:
             safeStack[k] = stackList   # Only assign if list NOT empty!


    return {'dataDims':dims, 'dataDimsUS':dimsUS, 'refDims':refDims, 'refDimsUS':refDimsUS,
            'shared':sharedDims, 'extra':extraDims, 'extraUS':extraDimsUS,
            'invalid':invalidDims, 'invalidUS':invalidDimsUS,
            'stacked':stackedDims, 'stackedMap':stackedDimsMap,
            'stackedShared':sharedDimsStacked, 'stackedExtra':extraDimsStacked, 'stackedInvalid':invalidDimsStacked,
            'missing':missingDims, 'safeStack':safeStack}



# Restack Xarray from refDims
# def restack(data, refDims):
#     """
#     Restack Xarray to conform to refDims.
#
#     Wraps checkDims() and data.stack() for "safe" restacking even with missing dims.
#
#     """
#
#     dimsDict = checkDims(data,refDims)
#
#     data = data.stack(dimsDict['safeStack'])
#
#     return data

# Restack Xarray from refDims - IN PROGRESS!!!
def restack(data, refDims = None, conformDims = False,
             forceUnstack = True, addMissing = False, unstackExtra = False,
             dimMap = None, strDims = ['Cont','Targ','Total','Type'],
             verbose = True):
    """
    Restack Xarray to conform to refDims.

    Wraps checkDims() and data.stack() for "safe" restacking even with missing dims.

    Parameters
    ----------

    data : Xarray
        Data to restack.

    refDims : optional, str, list or dict. Default = None.
        Reference dimensions
        - If None, use self.attrs['dataType']
        - If string, use epsproc.util.listFuncs.dataTypesList()[refDims]['def'](sType='sDict')
        - If list, treat as unstacked ref dims (will do dim check and missing dims only).
        - If dict, treat as stacked dim mappings.

    conformDims : bool, default = False
        If True, conform stacked dims to ref as closely as possible, including adding missing dims and unstacking extra dims.
        This sets forceUnstack = True, addMissing = True and unstackExtra = True.
        Note that extra dims will not be deleted.

    forceUnstack : bool, default = True
        Unstack input DataArray before further manipulation if True.

    addMissing : bool, default = False
        Add missing (unstacked) dims if True.
        Note these are added as string or numerical coords, as defined by strDims.
        (Setting incorrect coord type can affect some selection functions, specifically lmplot())

    unstackExtra : bool, default = False
        Unstack any extra stacked dims.
        Note that the dims are not removed, just unstacked.

    dimMap : dict, optional, default = None
        NOT YET IMPLEMENTED
        Map for dim renaming.

    strDims : list, default = ['Cont','Targ','Total','Type']
        Settings for addMissing for string coord dims.

    verbose : bool, default = True
        Show extra output if True.


    Returns
    -------
    Xarray : with restacked dims.

    Dict : output from checkDims() used to define the restacking.


    TODO:

    - Fix coord type issues, maybe define in dataTypesList?
    - Logging for before & after dims?

    """

    daTest = data.copy()  # Might be overkill

    # Set options to conform dims
    if conformDims:
        forceUnstack = True
        addMissing = True
        unstackExtra = True

    # if refDims is None:
    #     try:
    #         refDims = dataTypesList()[daTest.attrs['dataType']]['def'](sType='sDict')
    #     except:
    #         print("Input Xarray missing self.attrs['dataType']. Please set or pass refDims as string, list or dict to fix.")
    #
    # if isinstance(refDims, str):
    #     refDims = dataTypesList()[refDims]['def'](sType='sDict')

    # Above now in function
    refDims = getRefDims(daTest, refType = refDims, sType='sDict')

    # dimsDict = checkDims(daTest,refDims)

    if forceUnstack:    # and dimsDict:  # TODO - check existing dims first, and only proceed if restacking required?
        daTest = daTest.unstack()  # May not be necessary? But will allow for more complex cases.

    dimsDict = checkDims(daTest,refDims)

    # Add missing dims (uses UNSTACKED dim ref)
    if addMissing:
        for dim in dimsDict['invalidUS']:
            # Set for int vs. str labelled dims.
            # This works with existing lmplot routine OK
            if dim in strDims:
                daTest = daTest.expand_dims({dim: ['U']})
            else:
                daTest = daTest.expand_dims({dim: [0]})
        #         daTest = daTest.expand_dims({dim: [0]})

            if verbose:
                print(f'Added dim {dim}')

        dimsDict = checkDims(daTest,refDims)  # Update dimsDict in case new dims are stacked

    # Unstack extra dims?
    if dimsDict['stackedExtra'] and unstackExtra:
        daTest = daTest.unstack(dimsDict['stackedExtra'])

    # Safe stack to match ref.
    if dimsDict['safeStack']:
        daTest = daTest.stack(dimsDict['safeStack'])

    return daTest, dimsDict


# Subselect from sharedDims
def subselectDims(data, refDims = []):
    """
    Subselect dims from shared dim dict.
    Check dimensions for a data array (Xarray) vs. a reference list.

    Used to set safe selection criteria in matEleSelector.
    """

    # Check dims
    dimSets = checkDims(data, refDims)

    # Subselect
    if isinstance(refDims,dict):
        # Return dim with only subselected keys, i.e. existing dims.
        return {k:v for k,v in refDims.items() if k in dimSets['shared']}

    else:
        if not isinstance(refDims, list):
            refDims = [refDims]  # Force to list.

        # return dimSets['shared']  # Return shared dim list only.
        return [item for item in refDims if item in dimSets['shared']]   # Return shared items with forced ordering to match inputs.
