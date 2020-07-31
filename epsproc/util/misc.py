# -*- coding: utf-8 -*-
"""
ePSproc convenience functions

Collection of small functions for sorting etc.


"""

# import numpy as np
import re
import itertools
import os
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


# Sort & group filenames
def fileListSort(fList, groupByPrefix=True, verbose=True):
    """
    Sort a list of file names, and group by prefix.

    Note: this currently assumes a file name schema whereby split('_')[0] picks the grouping string.
    """

    if natsortFlag:
        fListSorted = natsorted(fList)
    else:
        fListSorted = sorted(fList)

    prefixStr = ''
    if groupByPrefix:
        prefixStr = os.path.commonprefix(fListSorted)  # Find common prefix

        # Solution with itertools groupby
        # Adapted from https://stackoverflow.com/a/13368753
        groupedList = [list(v) for k,v in itertools.groupby(fListSorted,key=lambda x:x.replace(prefixStr,'').split('_')[0])]


    if verbose:
        print(prefixStr)

    return fListSorted, groupedList, prefixStr
