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
# import scipy.constants
#
# # Package fns.
# from epsproc.basicPlotters import molPlot

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

    Returns
    -------

    dictionary
        Containing:

        - stacked and unstacked dims
        - stacked dim mappings
        - intersection and differences vs. refDims


    TODO: check and order dims by size? Otherwise set return is alphebetical

    26/08/21 Added additional tests for stacked dims vs. ref.
             Bit messy, but left as-is to avoid breaking old code.
             In future should amalgamate stacked & unstacked tests & tidy output.

    23/08/21 Added stacked dim mapping output.
    11/05/21 Added handling for stacked dims.

    """

    # Force to list to avoid breaking *unpacking later
    if not isinstance(refDims, (list, dict)):
        refDims = [refDims]

    dims = data.dims # Set dim list - this excludes stacked dims
    dimsUS = data.unstack().dims  # Set unstaked (full) dim list

    stackedDims = list(set(dims) - set(dimsUS))

    stackedDimsMap = {k:data.indexes[k].names for k in stackedDims}  # Get stacked dim mapping from indexes (Xarray/Pandas)

    # Check ref vs. full dim list (unstacked dims)
    sharedDims = list(set(dimsUS)&{*refDims})  # Intersection
    extraDims = list(set(dimsUS) - {*refDims})  # Difference
    invalidDims = list({*refDims} - set(dimsUS))  # Missing

    # 26/08/21 - added additional tests for stacked dims vs. ref, but may want to amalgamate with above?
    # TODO: tidy output too - separate dicts for stacked and unstacked dims?
    # Check ref vs. full dim list (stacked dims), note also unstacked dims subtracted to avoid duplicated dims.
    sharedDimsStacked = list(set(dims)&{*refDims} - set(dimsUS))  # Intersection
    extraDimsStacked = list(set(dims) - {*refDims} - set(dimsUS))  # Difference
    invalidDimsStacked = list({*refDims} - set(dims) - set(dimsUS))  # Missing

    # Test also missing dims overally
    missingDims = list({*refDims} - set(dimsUS) - set(dims))  # Missing

    return {'dataDims':dims, 'dataDimsUS':dimsUS, 'refDims':refDims,
            'shared':sharedDims, 'extra':extraDims, 'invalid':invalidDims,
            'stacked':stackedDims, 'stackedMap':stackedDimsMap,
            'stackedShared':sharedDimsStacked, 'stackedExtra':extraDimsStacked, 'stackedInvalid':invalidDimsStacked,
            'missing':missingDims}


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
        return dimsSets['shared']  # Return shared dim list only.
