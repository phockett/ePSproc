# -*- coding: utf-8 -*-
"""
ePSproc IO functions.
=====================

Module for file IO and data parsing.

Main function: :py:func:`epsproc.IO.readMatEle`:

    readMatEle(fileIn = None, fileBase = None, fType = '.out'):

    Read ePS file(s) and return results as Xarray data structures containing matrix elements.
    File endings specified by fType, default .out.


History
-------
07/06/22        Added various improvements to writeXarray and readXarray functionality.
                See https://github.com/phockett/ePSproc/issues/8 for ongoing notes.

13/10/20        Adapted main function readMatEle() to use grouped lists for multi-file jobs, should be back-compatible if stackE = False set.

06/11/19        Added jobInfo and molInfo data structures, from ePS file via :py:func:`epsproc.IO.headerFileParse()` and :py:func:`epsproc.IO.molInfoParse()`.
                Still needs a bit of work, and may want to implement other (comp chem) libraries here.

14/10/19        Added/debugged read functions for CrossSecion segments.

27/09/19        Added read functions for EDCS segments.

17/09/19        Added read/write to/from netCDF files for Xarrays.
                Use built-in methods, with work-arounds for complex number format issues.

29/08/19        Updating docs to rst.

26/08/19        Added parsing for E, sym parameters from head of ePS file.
                Added error checking by comparing read mat elements to expected list.
                Changed & fixed Xarray indexing - matrix elements now output with dims (LM, Eke, Sym, mu, it, Type)
                Current code rather ugly however.

19/08/19        Add functions for reading wavefunction files (3D data)

07/08/19        Naming convention tweaks, and some changes to comments, following basic tests with Sphinx.

05/08/19    v1  Initial python version.
                Working, but little error checking as yet. Needs some tidying.


To do
-----
* Add IO for other file segments (only DumpIdy supported so far).
* Better logic & flexibility for file scanning.
* Restructure as class for brevity...?
* More sophisticated methods/data structures for job & molecule info handling.



"""

# Imports
import os
import re
import numpy as np
import pandas as pd
from io import StringIO
import xarray as xr
from datetime import datetime as dt # Import datetime.datetime for now() function

try:
    from pyevtk.hl import gridToVTK
except ImportError as e:
    if e.msg != "No module named 'pyevtk'":
        raise
    print('* pyevtk not found, VTK export not available. ')

# Package fns.
from epsproc.util import matEleSelector, dataGroupSel, matEdimList, BLMdimList, stringRepMap, conv_ev_atm, orb3DCoordConv
from epsproc.util.misc import fileListSort, restack
from epsproc.util.env import isnotebook

# Set HTML output style for Xarray in notebooks (optional), may also depend on version of Jupyter notebook or lab, or Xr
# See http://xarray.pydata.org/en/stable/generated/xarray.set_options.html
if isnotebook():
    # Try/except here to allow for different XR versions.
    try:
        xr.set_options(display_style = 'html')
    except ValueError:
        pass

# ***** Ancillary functions

# File parsing function - scan file for keywords & read segments
#   Following above idiomatic solution, with full IO
#   https://stackoverflow.com/questions/3961265/get-line-number-of-certain-phrase-in-file-python
def fileParse(fileName, startPhrase = None, endPhrase = None, comment = None, verbose = 0):
    """
    Parse a file, return segment(s) from startPhrase:endPhase, excluding comments.

    Parameters
    ----------
    fileName : str
        File to read (file in working dir, or full path)
    startPhrase : str, optional
        Phrase denoting start of section to read. Default = None
    endPhase : str or list, optional
        Phrase denoting end of section to read. Default = None
    comment : str, optional
        Phrase denoting comment lines, which are skipped. Default = None
    verbose : int, optional, default = 1
        Level of verbosity in output.
        - 0 no printed output
        - 1 print summary info only
        - 2 print detailed info

    Returns
    -------
    list
        [lineStart, lineStop], ints for line #s found from start and end phrases.
    list
        segments, list of lines read from file.

    All lists can contain multiple entries, if more than one segment matches the search criteria.

    """

    lineStart = []    # Create empty list to hold line #s
    lineStop = []     # Create empty list to hold line #s
    segments = [[]]   # Possible to create empty multi-dim array here without knowing # of segments? Otherwise might be easier to use np textscan functions
    readFlag = False
    n = 0

    # Force list to ensure endPhase is used correctly for single phase case (otherwise will test chars)
    if type(endPhrase) is str:
        endPhrase = [endPhrase]

    # Open file & scan each line.
    with open(fileName,'r') as f:
        for (i, line) in enumerate(f):  # Note enumerate() here gives lines with numbers, e.g. fullFile=enumerate(f) will read in file with numbers
            i = i + 1  # Offset for file line numbers (1 indexed)
            # If line matches startPhrase, print line & append to list.
            # Note use of lstrip to skip any leading whitespace.
            # if startPhrase in line:
            if line.lstrip().startswith(startPhrase):
                if verbose > 1:
                    print('Found "', startPhrase, '" at line: ', i)

                lineStart.append(i)

                readFlag = True

            # Read lines into segment[] until endPhrase found
            if readFlag:
                # Check for end of segment (start of next Command sequence)
                if endPhrase and ([line.startswith(endP) for endP in endPhrase].count(True) > 0):  # This allows for multiple endPhases
                                                                                                    # NOTE: this will iterate over all chars in a phrase if a single str is passed.
                    # Log stop line and list
                    lineStop.append(i)
                    readFlag = False

                    # Log segment and create next
                    segments.append([])
                    n += 1

                    continue            # Continue will skip rest of loop

                 # Check for comments, skip line but keep reading
                elif comment and line.startswith(comment):
                    continue            # Continue will skip rest of loop

                segments[n].append([n, i, line])    # Store line if part  of defined segment

    if verbose:
        # print('Found {0} segments.'.format(n+1))
        print(f'*** IO.fileParse() found {n} segments with \n\tStart: {startPhrase}\n\tEnd: {endPhrase}.')

    return ([lineStart, lineStop], segments) # [:-1])


# Simple wrapper for general fileParse function, ePS dumpIdy segments
def dumpIdyFileParse(fileName, verbose = 1):
    """
    Parse an ePS file for dumpIdy segments.

    Parameters
    ----------
    fileName : str
        File to read (file in working dir, or full path)

    verbose : bool or int, optional
        If true, print segment info.

    Returns
    -------
    list
        [lineStart, lineStop], ints for line #s found from start and end phrases.
    list
        dumpSegs, list of lines read from file.

    Lists contain entries for each dumpIdy segment found in the file.

    """

    startPhrase = "DumpIdy - dump"
    endPhrase = ["+ Command", "Time Now"]  # In this case may have multiple end phrases

    (lines, dumpSegs) = fileParse(fileName, startPhrase, endPhrase, verbose = verbose) # , '>')

    # NOTE - with current code dumpSegs has one final blank segment
    if verbose:
        print('Found {0} dumpIdy segments (sets of matrix elements).'.format(len(dumpSegs) - 1))

    return lines, dumpSegs[:-1]

# Simple wrapper for general fileParse function, ePS EDCS segments
def EDCSFileParse(fileName, verbose = 1):
    """
    Parse an ePS file for EDCS segments.

    Parameters
    ----------
    fileName : str
        File to read (file in working dir, or full path)

    verbose : bool or int, optional
        If true, print segment info.

    Returns
    -------
    list
        [lineStart, lineStop], ints for line #s found from start and end phrases.
    list
        dumpSegs, list of lines read from file.

    Lists contain entries for each dumpIdy segment found in the file.

    """

    startPhrase = "EDCS - differential cross section program"
    endPhrase = ["+ Command", "Time Now"]  # In this case may have multiple end phrases

    (lines, dumpSegs) = fileParse(fileName, startPhrase, endPhrase, verbose = verbose) # , '>')

    # NOTE - with current code dumpSegs has one final blank segment
    if verbose:
        print('Found {0} EDCS segments (sets of scattering results).'.format(len(dumpSegs) - 1))

    return lines, dumpSegs[:-1]

# Simple wrapper for general fileParse function, ePS GetCro/CrossSection segments
def getCroFileParse(fileName, verbose = 1):
    """
    Parse an ePS file for GetCro/CrossSection segments.

    Parameters
    ----------
    fileName : str
        File to read (file in working dir, or full path)

    verbose : bool or int, optional
        If true, print segment info.

    Returns
    -------
    list
        [lineStart, lineStop], ints for line #s found from start and end phrases.
    list
        dumpSegs, list of lines read from file.

    Lists contain entries for each dumpIdy segment found in the file.

    """

    # startPhrase = "CrossSection - compute photoionization cross section"  # For full file segment
    startPhrase = "COMPOSITE CROSS SECTIONS AT ALL ENERGIES"  # Final table of XSect values
    endPhrase = ["+ Command", "Time Now"]  # In this case may have multiple end phrases

    (lines, dumpSegs) = fileParse(fileName, startPhrase, endPhrase, verbose = verbose) # , '>')

    # NOTE - with current code dumpSegs has one final blank segment
    if verbose:
        print('Found {0} CrossSection segments (sets of results).'.format(len(dumpSegs) - 1))

    return lines, dumpSegs[:-1]



# Simple wrapper for general fileParse function, check ScatEng and return Eke list
def scatEngFileParse(fileName, verbose = 1):
    """
    Parse an ePS file for ScatEng list.

    Parameters
    ----------
    fileName : str
        File to read (file in working dir, or full path)

    verbose : bool or int, optional
        If true, print segment info.

    Returns
    -------
    list
        ekeList, np array of energies set in the ePS file.

    Lists contain entries for each dumpIdy segment found in the file.

    """
    startPhrase = "ScatEng"
    endPhrase = "#"
    comment = "+"

    (lines, dumpSegs) = fileParse(fileName, startPhrase, endPhrase, comment, verbose = verbose) # , '>')

    # Grab E list, assuming just first segment scanned is relevant
    ekeList = np.genfromtxt(StringIO(dumpSegs[0][0][2][7:]))

    # print('Expecting {0} energy points.'.format(len(ekeList)))  # For np.array, this will fail for singleton dim.
    if verbose:
        print('Expecting {0} energy points.'.format(ekeList.size))

    return ekeList

# Simple wrapper for general fileParse function, check symmetries and return list
def symFileParse(fileName, verbose = 1):
    """
    Parse an ePS file for scattering symmetries.

    Parameters
    ----------
    fileName : str
        File to read (file in working dir, or full path)

    verbose : bool or int, optional
        If true, print segment info.

    Returns
    -------
    list
        symSegs, raw lines from the ePS file.

    Lists contain entries for each ScatSym setting found in file header (job input).

    """
    startPhrase = "ScatSym"
    endPhrase = ["FileName", "\n"]
    comment = "+"

    (lines, symSegs) = fileParse(fileName, startPhrase, endPhrase, comment, verbose = verbose) # , '>')

    # Grab E list, assuming just first segment scanned is relevant
    if verbose:
        print('Expecting {0} symmetries.'.format(len(symSegs) - 1))

    return symSegs[:-1]


# Parse digits from a line using re
# https://stackoverflow.com/questions/4289331/how-to-extract-numbers-from-a-string-in-python
def parseLineDigits(testLine):
    """
    Use regular expressions to extract digits from a string.
    https://stackoverflow.com/questions/4289331/how-to-extract-numbers-from-a-string-in-python

    """
    return re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?", testLine)



# ************** Header info parsing functions

# Simple wrapper for general fileParse function, extract ePS file header & input, and parse
def headerFileParse(fileName, verbose = True):
    """
    Parse an ePS file for header & input job info.

    Parameters
    ----------
    fileName : str
        File to read (file in working dir, or full path)

    verbose : bool, default True
        Print job info from file header if true.

    Returns
    -------
    jobInfo : dict
        Dictionary generated from job details.

    TO DO
    -----
    - Tidy up methods - maybe with parseDigits?
    - Tidy up dict output.

    """

    startPhrase = "ePolyScat Version" # Read from top of file
    endPhrase = ["+ End of input reached"]  # In this case only have a single end phrases, but need to pass as list to avoid iterating over phrase

    (lines, dumpSegs) = fileParse(fileName, startPhrase, endPhrase, verbose = verbose) # , '>')

    # NOTE - with current code dumpSegs has one final blank segment
    # print('Read {0} dumpIdy segments (sets of matrix elements).'.format(len(dumpSegs) - 1))

    # Parse info to dict - bit ugly, assumes fixed format for start lines
    # This might be workable if keys are set from ePS source/web list?
    # jobKeys = ['ePolyScat', 'Authors', 'http', 'cite', 'Starting', 'Using']
    # jobInfo = {}
    #
    # for line in dumpSegs[0]:
    #     for key in jobKeys:
    #         if line[2].startswith(key):
    #             jobInfo[key] = line[2].strip()

    # Generate from file directly - good for key:value pairs, but might mangle prose
    jobInfo = {}
    jobInfo['comments'] = []    # Keep comments in this list.

    # Loop over lines, split and sort if possible - UGLY!
    for n, line in enumerate(dumpSegs[0]):
    #     elements = dumpSegs[0][n][2].strip().split()
        elements = line[2].strip().split()

        # print(elements)
        if len(elements)>0:
            if elements[0].startswith('#'):
                jobInfo['comments'].append(line[2].strip())
            elif elements[0].startswith('Orb'):
                jobInfo[elements[0]] = np.asarray(parseLineDigits(dumpSegs[0][n+1][2].strip())).astype('int')  # Set next line for OrbOccs
            else:
                if len(elements) == 2:
                    jobInfo[elements[0]] = elements[1]  # Case for data record value assignments

                # Check for lines with comments in
                elif '#' in line[2]:
                    test = line[2].strip().split('#')  # Split at comment
                    testStart = test[0].split()
                    if len(testStart) == 2:
                        jobInfo[testStart[0]] = testStart[1]  # Case for data record assignments

                    # This might be redundant... but keep for safety.
                    else:
                        jobInfo[elements[0]] = line[2].strip().split('#')   # For other cases keep full line, split at comments

                # For all other cases keep full line, split at comments
                else:
                    jobInfo[elements[0]] = line[2].strip().split('#')

    # Print jobInfo
    if verbose:
        print('*** Job info from file header.\n')
        [print(line.strip('#')) for line in jobInfo['comments'][0:4]]

    return jobInfo

# Simple wrapper for general fileParse function, extract ePS molecular info, and parse
def molInfoParse(fileName, verbose = True):
    """
    Parse an ePS file for input molecule info.

    Parameters
    ----------
    fileName : str
        File to read (file in working dir, or full path)

    verbose : bool, default True
        Print job info from file header if true.

    Returns
    -------
    molInfo : dict
        Dictionary with atom & orbital details.

    Notes
    -----

    Only tested for Molden input (MoldenCnv2006).

    """

    # Extract mol segment
    startPhrase = "Selecting orbitals" # Read from top of file
    endPhrase = ["+ "]  # In this case only have a single end phrases, but need to pass as list to avoid iterating over phrase
    (lines, dumpSegs) = fileParse(fileName, startPhrase, endPhrase, verbose = verbose)

    # Basic parsing to lists
    orbList = []
    atomList = []

    for line in dumpSegs[0]:
        if line[2].startswith('Selecting'):
            orbList.append(line[2].strip().lstrip('Selecting '))
        if line[2].startswith('Z'):
            atomList.append(line[2].strip())

    orbList = orbList[1:]

    if verbose:
        print('\n*** Found orbitals')
        print(*orbList, sep='\n')
        print('\n*** Found atoms')
        print(*atomList, sep='\n')

    # Sort orbs to np.array
    orbTable = []
    [orbTable.append(parseLineDigits(orb)) for orb in orbList]
    orbTable = np.asarray(orbTable).astype('float')

    # # Convert to Xarray - now moved below to include additional orb info.
    # cols = ['N', 'EH', 'Occ', 'E']
    # orbTableX = xr.DataArray(orbTable[:,1:], coords = {'orb':orbTable[:,0], 'props':cols}, dims=['orb','props'])

    # 26/03/21 - moved labels here with basic dim check to allow for different output in ePS build 3885d87
    if 'SymOrb' in orbList[0]:
        orbListCols = ['N', 'SymOrb', 'EH', 'Occ', 'E']
        EHind = 3
    else:
        orbListCols = ['N', 'EH', 'Occ', 'E']
        EHind = 2

    orbTable = np.c_[orbTable, conv_ev_atm(orbTable[:,EHind], to='ev')]  # Convert to eV and add to array
    EVind = orbTable.shape[1] - 2  # Set E index for Xarray conv later.

    # Sort coords to np.array
    atomTable = []
    [atomTable.append(parseLineDigits(atom)) for atom in atomList]
    atomTable = np.asarray(atomTable).astype('float')

    # Convert to Xarray
    cols = ['Z', 'Zs', 'x', 'y', 'z']
    atomTableX = xr.DataArray(atomTable, coords = {'atom':np.arange(1,atomTable.shape[0]+1), 'props':cols}, dims=['atom','props'])


    #*** Additional info from other segments

    # Read SymProd output, this lists orbs by symmetry & degenerate groups
    startPhrase = "SymProd - Construct products of symmetry types" # Read from top of file
    endPhrase = ["Time Now ="]  # In this case only have a single end phrases, but need to pass as list to avoid iterating over phrase
    (lines, dumpSegs) = fileParse(fileName, startPhrase, endPhrase, verbose = verbose)

    # Convert lines to useful values
    orbSet = []
    orbList = []
    orbSymList = []

    # Parse for specific line beginnings only ... For grouping, easier to explicitly check values?

    # Parse 'Set' and 'Orbital' lines
    lineStart = ['Set', 'Orbital']
    for item in dumpSegs[0]:
        if(item[2].startswith(lineStart[0])):
            orbSet.append(np.asarray(parseLineDigits(item[2])).astype('int'))

        if(item[2].startswith(lineStart[1])):
            temp = item[2].split()

            # Assign values
            if len(temp) == 12:
                orbList.append([temp[4], temp[7], temp[11], orbSet[-1][0], orbSet[-1][1]])
                orbSymList.append(temp[10])

    # Convert final results to np.array
    orbList = np.asarray(orbList).astype('int')
    temp = np.c_[orbTable[:,1:], orbList[:,1:]]

    # Parse ExpOrb output to check orbital expansion convergence
    startPhrase = "ExpOrb - Single Center Expansion Program" # Read from top of file
    endPhrase = ["Time Now ="]  # In this case only have a single end phrases, but need to pass as list to avoid iterating over phrase
    (lines, dumpSegs) = fileParse(fileName, startPhrase, endPhrase, verbose = verbose)

    # Parse 'Orbital' lines
    # Note - these are orbital groups
    orbIntList = []
    for item in dumpSegs[0]:
        if(item[2].startswith('Orbital')):
            # temp = item[2].split()
            # orbIntList.append(parseLineDigits(item[2]))
            orbMetrics = parseLineDigits(item[2])
            orbIntList.append([orbMetrics[0], orbMetrics[-1]])  # Set explicitly here, length may change depending on whether symmetries have digits in!

    orbIntList = np.asarray(orbIntList).astype('float')

    # Assign to main table by orbGrp
    # - THIS IS HORRIBLE, should do with Xarrays directly, or just better sorting code, but not happening right now. It's Monday.
    ogInt = np.zeros(temp.shape[0])
    for n, row in enumerate(temp):
        ind = row[6] == orbIntList[:,0]  # Find corresponding row
        ogInt[n]=orbIntList[ind, 1]  # Set int value

    temp = np.c_[temp, ogInt]

    #*** Set as Xarray
    # TODO: make this neater/better!
    # cols = ['N', 'EH', 'Occ', 'E', 'Type', 'NOrbGrp', 'OrbGrp', 'GrpDegen', 'NormInt']
    orbListCols.extend(['Type', 'NOrbGrp', 'OrbGrp', 'GrpDegen', 'NormInt'])  # 26/03/21 adapted for new ePS build
    orbTableX = xr.DataArray(temp, coords = {'orb':orbTable[:,0].astype(int), 'props':orbListCols}, dims=['orb','props'])  # TODO: fix type here - changes to float in Xarray...?
    orbTableX['Sym'] = ('orb', orbSymList)  # Add symmetries, E and OrbGrp as non-dim coords - possibly not a great thing to do...?  Not easy to index with these
    orbTableX['E'] = ('orb', temp[:,EVind])
    orbTableX['OrbGrp'] = ('orb', orbList[:,3])

    # Alternatively - add as dim coords. This duplicates data in this case, would need to sort by line or reshape array to get this working (cf. matrix elements)
    # orbTableX = orbTableX.expand_dims({'Sym':orbSymList}) # , 'Eke':[attribs[0][1]]})

    # Compile to dict
    molInfo = {'atomList':atomList, 'atomTable':atomTableX, 'orbList':orbList, 'orbTable':orbTableX}

    return molInfo


# ************* DumpIdy parsing

# Parse a DumpIdy segment
#TODO: More attribs, and convert attribs to dict, or other more structured format.
#TODO: Error checking for missing matrix elements or IO issues - see Matlab code.
#TODO: Convert attribs to dict for simpler assignments later
def dumpIdySegParse(dumpSeg):
    """
    Extract values from dumpIdy file segments.

    Parameters
    ----------
    dumpSeg : list
        One dumpIdy segment, from dumpSegs[], as returned by :py:func:`epsproc.IO.dumpIdyFileParse()`

    Returns
    -------
    np.array
        rawIdy, array of matrix elements, [m,l,mu,ip,it,Re,Im]
    list
        attribs, list [Label, value, units]

    Notes
    -----
    Currently this is a bit messy, and relies on fixed DumpIDY format.
    No error checking as yet.
    Not yet reading all attribs.

    Example
    -------

    >>> matEle, attribs = dumpIdySegParse(dumpSegs[0])

    """

    # Use lists to collect data, and convert format at the end
    attribs = []
    rawIdy = []

    # print(len(dumpSeg))
    # Parse data block
    # Use native python, or np.genfromtxt
    for testLine in dumpSeg[12:-1]:

#        vals = parseLineDigits(testLine[2])
#        lineOut = []
#        [lineOut.append(float(val)) for val in vals] # Convert to float
#        rawIdy.append(lineOut)

        tmp=np.genfromtxt(StringIO(testLine[2]))
        rawIdy.append(tmp)

    # Parse header lines (structured, selected by line #)
    attribs.append(['E', np.float(parseLineDigits(dumpSeg[3][2])[0]), 'eV'])
    attribs.append(['Ehv', np.float(parseLineDigits(dumpSeg[4][2])[0]), 'eV'])
    SF = np.genfromtxt(parseLineDigits(dumpSeg[5][2]))
    SF = SF[0] + SF[1]*1j
    attribs.append(['SF', SF, 'sqrt(MB)'])
    attribs.append(['Lmax', np.int(parseLineDigits(dumpSeg[10][2])[0]), ''])

    # Parse symmetries - multiple .spilt(), or use re?
    # attribs.append(['Syms', parseLineDigits(dumpSeg[7][2]), ''])
    symList = dumpSeg[7][2].split()
    # Append as nested lists
#    sym = []
#    [sym.append([symList[n], symList[n+1].strip('=')]) for n in range(1,6,2)]
#    attribs.append(['Syms', sym, ''])
    # Append as indendent lists
    [attribs.append([symList[n], symList[n+1].strip('=')]) for n in range(1,6,2)]

    attribs.append(['QNs', dumpSeg[11][2].split(), ''])

    return np.asarray(rawIdy), attribs


# Functional form for parsing full set of mat elements and putting in xarray
def dumpIdySegsParseX(dumpSegs, ekeListUn, symSegs, verbose = 1):
    """
    Extract data from ePS dumpIdy segments into usable form.

    Parameters
    ----------
    dumpSegs : list
        Set of dumpIdy segments, i.e. dumpSegs, as returned by :py:func:`epsproc.IO.dumpIdyFileParse()`

    ekeListUn : list
        List of energies, used for error-checking and Xarray rearraging, as returned by :py:func:`epsproc.IO.scatEngFileParse()`

    verbose : bool, default True
        Print job info from file header if true.

    Returns
    -------
    xr.array
        Xarray data array, containing matrix elements etc.
        Dimensions (LM, Eke, Sym, mu, it, Type)

    int
        Number of blank segments found.


    Example
    -------

    >>> data = dumpIdySegsParseX(dumpSegs)

    """

    # Logging
    dataList = []
    ekeList = []
    blankSegs = 0

    blankMatE = 0
    blankSymList = []

    # Loop over DumpIdy segments, extract data & reformat
    # If blank, skip parser and append blankSegs.
    for dumpSeg in dumpSegs:
        if len(dumpSeg)>6:
            segBlock, attribs = dumpIdySegParse(dumpSeg)

            # Test also on returned sebBlock, this can be empty in some cases (matrix elements below threhold, which creates stub in output file but contains no values)
            # if not segBlock.any():
                # segBlock = np.zeros([1,7])  # With zeros - gives index issues later
                # segBlock = np.full([1,7], fill_value = np.nan)  # With NaNs - gives issues later however.

            if segBlock.any():
                dataList.append([segBlock[:,0:5].T, segBlock[:,5]+1j*segBlock[:,6], attribs])
                ekeList.append(attribs[0][1])

            # Skip and treat as blankSeg in this case, as per default case below
            # HOWEVER - may want to propagate missing symmetry case here?
            # ACTUALLY - propagate as blankMatE, and skip Eke listing.
            else:
                blankMatE += 1
                # ekeList.append(attribs[0][1])
                blankSymList.append(attribs[4][1])

                # ekeList.append(np.nan)  # This adds multiple np.nan entries which propagate through np.unique() later.

            # else:
            #     # dataList.append([*(6*[np.nan]), attribs])  # Set list of nans for empty case.
            #     dataList.append([*(6*[np.nan]), attribs])


            # Switch l,m - with advanced indexing, other methods faster...? https://stackoverflow.com/questions/4857927/swapping-columns-in-a-numpy-array
            # Now done later via pd.MultiIndex.swaplevel()
            # dataList[-1][0][[0,1],:] = dataList[-1][0][[1,0],:]

            # dataList.append([segBlock[:,0:5].T, segBlock[:,5]+1j*segBlock[:,6], attribs])
            # dataList.append([segBlock[:,0:5], segBlock[:,5]+1j*segBlock[:,6], attribs])



        else:
            blankSegs += 1
            ekeList.append(np.nan)

    # Check energies vs. input list, and number of symmetries
    ekeTest = np.unique(ekeList)
    if ekeTest.size != ekeListUn.size:
        print("*** Warning: Found {0} energies, expected {1}".format(ekeTest.size,ekeListUn.size))

    # Check here according to expected input, but should also be logged in blankSegs above.
    if len(ekeList) != ekeTest.size * len(symSegs):
        print("*** Warning: Missing records, expected {0}, found {1}.".format(ekeTest.size * len(symSegs),len(ekeList)))

    if blankMatE > 0:
        print("*** Warning: Found {0} blank sets of matrix elements, symmetries {1}".format(blankMatE, np.unique(blankSymList)))

    #**** Convert to xarray - ugly loop version, probably a better way to do this!
    #TODO Should:
    #   - be able to loop more cleanly over attribs - set as dict?
    #   - integrate with above loop
    #   - compile into dataSet or dataArray more directly...?
    #   - Check and compile against Eke list (muliple symetries per Eke), and keep this as a separate coord.
    dataArrays = []
    dataSym = []
    ekeVal = ekeList[0]

    if verbose:
        print('\nProcessing segments to Xarrays...')

    eLoop = 0   # For logging E vs. sym segments and looping output.
    for n, data in enumerate(dataList):
        attribs = data[2]

        # V1 - EASIEST WAY, but leads to dimensional issues later!
        #TODO: consider setting mu as a separate dim. Maybe also (ip,it)...?
        # QNs = pd.MultiIndex.from_arrays(data[0].astype('int8'), names = attribs[-1][1][0:-1])
        # QNs = QNs.swaplevel(0, 1)  # Switch l,m indexes
        # Esyms = pd.MultiIndex.from_arrays([np.array(attribs[0][1]), [attribs[3][1], attribs[5][1]]], names=['E', 'Sym'])
        # pd.MultiIndex.from_tuples([(np.array(attribs[0][1]), [attribs[3][1], attribs[5][1]])], names=['E', 'Sym'])
        # Esyms = pd.MultiIndex.from_tuples([(attribs[0][1],attribs[4][1],attribs[5][1],attribs[6][1])],names=[attribs[0][0],attribs[4][0],attribs[5][0],attribs[6][0]])
        #dataArrays.append(xr.DataArray(data[1], coords={'ES': Esyms, 'LM':QNs}, dims = ['ES','LM']))

        # AH - issue is number of labels - can't lable singleton dim it seems, but can expand
        #TODO: consider setting E as a separate dim, will be singleton for each set of syms. Might make more sense for later manipulations (sum over sym or E).
        # tmp = xr.DataArray(np.asarray(data[1]), coords={'LM':QNs}, dims = ['LM'])
        # tmp = tmp.expand_dims({'Sym':Syms, 'Eke':[attribs[0][1]]})  # This is OK, but still ties Eke and Sym coords (same number of elements)
        # # tmp = tmp.expand_dims({'Sym':Syms})


        # # v2 - separate dims, (LM, E, Syms, ip, it)
#        LM = pd.MultiIndex.from_arrays(data[0][0:2,:].astype('int8'), names = attribs[-1][1][0:2])
#        LM = LM.swaplevel(0, 1)  # Switch l,m indexes
#        mu = data[0][2,:]
#        ip = data[0][3,:]
#        it = data[0][4,:]
#        Syms = pd.MultiIndex.from_tuples([(attribs[4][1],attribs[5][1],attribs[6][1])],names=[attribs[4][0],attribs[5][0],attribs[6][0]])
#
        # tmp =  xr.DataArray(np.asarray(data[1]), coords={'LM':LM, 'mu':mu, 'ip':ip, 'it':it}, dims = ['LM','mu','it','ip'])
        # This works... but still keeps full lenght for additional label coords.
        # Should be able to reduce on these... but can't work out how (tried groupby etc.)
#        tmp =  xr.DataArray(np.asarray(data[1]), coords={'LM':LM}, dims = ['LM'])
#        tmp = tmp.expand_dims({'Sym':Syms, 'Eke':[attribs[0][1]]})
#        tmp = tmp.expand_dims({'mu':mu, 'ip':ip, 'it':it})

        # *** v2b - assign as v1/v2, then sort Xarrays before restacking
        # This works, but note assumption of stacking order (E, then Syms)
        QNs = pd.MultiIndex.from_arrays(data[0].astype('int8'), names = attribs[-1][1][0:-1])
        QNs = QNs.swaplevel(0, 1)  # Switch l,m indexes
        Syms = pd.MultiIndex.from_tuples([(attribs[4][1],attribs[5][1],attribs[6][1])],names=[attribs[4][0],attribs[5][0],attribs[6][0]])

        # Original code - set according to LM, then expand dims.
        tmp = xr.DataArray(np.asarray(data[1]), coords={'LM':QNs}, dims = ['LM'])
        tmp = tmp.expand_dims({'Sym':Syms, 'Eke':[attribs[0][1]]})

        # New code 24/10/19 - set all coords, including non-dim coords such as SF and Ehv
        # Should be able to set in one call... but get dim issues here, probably with multiple singletons.
        # tmp = xr.DataArray(np.asarray(data[1]), coords={'LM':QNs, 'Sym':Syms, 'Eke':[attribs[0][1]], 'Ehv':[attribs[1][1]], 'SF':[attribs[2][1]]}, dims = ['LM', 'Sym', 'Eke'])
        # Setting as per previous code, then adding addtional singleton non-dimensional coords seems to be OK however...
        tmp['Ehv']=(('Eke'),[attribs[1][1]])  # Link to single dim... should be OK?
        # tmp['SF']=(('Eke'),[attribs[2][1]])  # Link to single dim... should be OK?
        tmp['SF']=(('Eke','Sym'),np.array(attribs[2][1]).reshape(1,1))  # Link to Eke & sym (i.e. single set of matE)

        # Assign any other attributes - note that some attributes may be dropped when combining arrays below
        for a in attribs:
            tmp.attrs[a[0]] = a[1] # Currently set without units, multiple values here give combine issues below.

        # dataArrays.append(tmp)

#        tmp = matEleGroupDimX(tmp)  # Broken?
        dataArrays.append(matEleGroupDimX(tmp))
        #TODO: Fix scale-factor propagation... this is currently dropped as an inconsistent attrib value.
        # Something like: dataSF.append(dataArrays[-1].SF)
        # dataArrays.append(matEleGroupDimXnested(tmp.copy()))  # Broken...?
        # dataArrays.append(tmp)

        #TODO: UGLY - need to check and combine according to number of Syms (check matlab code)
        # ALSO NOT WORKING PROPERLY - might be issue with equality?
        # BETTER TO JUST SET 2D arrays here!
#        if n == (len(ekeList)/2 - 1):
#            dataSym.append(xr.combine_nested(dataArrays, concat_dim=['Eke']))
#            dataArrays = []
#
#        if n == len(ekeList)-1:
#            dataSym.append(xr.combine_nested(dataArrays, concat_dim=['Eke']))
#
#        if n == (len(ekeList)/2 - 1):
#            # dataSym.append(xr.combine_nested(dataArrays, concat_dim=['Eke']))
#            dataArrays1 = dataArrays.copy()
#            dataArrays = []

        # if n == (len(ekeList)/2 - 1):
        #     dataSym.append(dataArrays)
        #     dataArrays = []
        #
        # if n == len(ekeList)-1:
        #     dataSym.append(dataArrays)

        # if (n > 1) and (attribs[0][1] == ekeVal):
        #     print('Found {0} energies'.format(n+1))
        #     dataSym.append(xr.combine_nested(dataArrays, concat_dim=['Eke']))


        # Now with loop on known energy list.
        if eLoop == ekeListUn.size-1:
            dataSym.append(xr.combine_nested(dataArrays, concat_dim=['Eke']))
            dataArrays = []
            eLoop = 0
        else:
            eLoop = eLoop+1


        # **** v3 manually sort data first...
        # SEE CODE BELOW, matEleGroupDim()

        # Assign any other attributes - note that some attributes may be dropped when combining arrays below
        # for a in attribs:
        #     tmp.attrs[a[0]] = a[1] # Currently set without units, multiple values here give combine issues below.
        #
        # dataArrays.append(tmp)

        # Stack by syms (per eke) as necessary
#        if ekeList[n] == ekeVal:
#            # daOut = xr.combine_nested(dataArrays, concat_dim=['Sym'])
#            dataArrays.append(tmp)
#        else:
#            dataSym.append(xr.combine_nested(dataArrays, concat_dim=['Sym']))
#            dataArrays = []
#            dataArrays.append(tmp)
#            ekeVal = ekeList[n]

    # Combine to single xarray
    # Note xarray > v0.12.1
    # daOut = xr.combine_nested(dataSym, concat_dim=['Sym'])
    daOut = xr.combine_nested(dataSym, concat_dim=['Sym'])
    # daOut = xr.combine_nested([dataArrays1, dataArrays], concat_dim=['Sym','Eke'])
    # daOut = xr.combine_nested(np.array(dataArrays), concat_dim=['Sym','Eke'])  # Dim issues here - try np.array... NOPE.
    # daOut = xr.combine_nested(dataArrays, concat_dim=[None])
    # daOut = xr.merge(dataArrays)
    # daOut = xr.combine_by_coords(dataArrays)
    # daOut = dataArrays
    # daOut = daOut.expand_dims({'Eke':[attribs[0][1]]})

#    # v3 - Sort data before putting into Xarray
#    # NOW REPLACED ABOVE by sorting of Xarrays - code here may be faster, but less robust.
#    dataArrays = []
#    for data in dataList:
#        attribs = data[2]
#        tmp = matEleGroupDim(data)
#        dataArrays.append(tmp)
#
#    # Can only concat over it at the moment due to duplicate values issue
#    # THIS WILL BREAK LATER!!!!
#    # TODO: fix use of Xarrays and dimension issues.
#    daOut = xr.combine_nested(dataArrays, concat_dim=['it'])

    # Set any other global attribs
    daOut.attrs['dataType'] = 'matE'    # Set dataType for use later.

    # SF testing vs. symmetry - remove Sym dim if not necessary, for easier computation later!
    # Added 25/10/19, not yet well tested.
    if daOut.SF.diff('Sym').max().pipe(np.abs) < np.finfo(complex).eps:
        daOut['SF'] = ('Eke', daOut.SF.values[:,0])


    return daOut.transpose(), blankSegs     # NOTE transpose to reverse ordering of dims.



# ************* EDCS parsing

# Parse a EDCS segment (roughly)
# def EDCSSegParse(dumpSeg):
#     """
#     Extract values from EDCS file segments.
#
#     Parameters
#     ----------
#     dumpSeg : list
#         One EDCS segment, from dumpSegs[], as returned by :py:func:`epsproc.IO.EDCSFileParse()`
#
#     Returns
#     -------
#     np.array
#         EDCS, array of scattering XS, [theta, Cross Section (Angstrom^2)]
#     list
#         attribs, list [Label, value, units]
#
#     Notes
#     -----
#     Currently this is a bit messy, and relies on fixed EDCS format.
#     No error checking as yet.
#     Not yet reading all attribs.
#
#     Example
#     -------
#
#     >>> EDCS, attribs = EDCSSegParse(dumpSegs[0])
#
#     """
#     # Use lists to collect data, and convert format at the end
#     attribs = []
#     EDCS = []
#
#     attribs.append(['E', np.float(parseLineDigits(dumpSeg[13][2])[0]), 'eV'])
#
#     # For each line convert to float - bit ugly, but works
#     for testLine in dumpSeg[67:]:
#         EDCS.append(np.genfromtxt(StringIO(testLine[2])))
#
#     return np.asarray(EDCS), attribs


# Updated function 07/03/22 - additional line checks to allow for variable-length output & debug output with verbose
def EDCSSegParse(dumpSeg, verbose = False):
    """
    Extract values from EDCS file segments.

    Parameters
    ----------
    dumpSeg : list
        One EDCS segment, from dumpSegs[], as returned by :py:func:`epsproc.IO.EDCSFileParse()`

    verbose : bool, int, default = False
        Print additional info during run.

    Returns
    -------
    np.array
        EDCS, array of scattering XS, [theta, Cross Section (Angstrom^2)]
    list
        attribs, list [Label, value, units]

    Notes
    -----

    Only reads (theta,I) information from an "EDCS - differential cross section program" segment.
    Currently this is a bit messy, and relies on fixed EDCS format.
    No error checking as yet.
    Not yet reading all attribs.

    Example
    -------

    >>> EDCS, attribs = EDCSSegParse(dumpSegs[0])

    """

    # Use lists to collect data, and convert format at the end
    attribs = []
    EDCS = []

    DCSFlag = False  # Use this to turn on output for DCS section of segment only.

    for testLine in dumpSeg:
    #     print(line)

        # Get energy, read line "Energy to compute the EDCS at (eV) =      X"
        # Note there may also be another energy line with other units provided.
        if testLine[-1].strip().startswith('Energy to compute'):
            if verbose:
                print(f"***E line: {testLine}")

            attribs.append(['E', np.float(parseLineDigits(testLine[2])[0]), 'eV'])

        # Use "Ang" or "Differential Cross Section" as keywords to flag beginning of (theta,I) data
        if testLine[-1].strip().startswith('Differential Cross Section'):
            if verbose:
                print(f"*** Seg {testLine[0]}, line {testLine[1]} DCS start")

            DCSFlag = True


        if testLine[-1].strip().startswith("Time Now"):
            if verbose:
                print(f"Seg {testLine[0]}, line {testLine[1]} DCS end")

            DCSFlag = False

        if DCSFlag:
            # For each line convert to np - bit ugly, but works
            testVals = np.genfromtxt(StringIO(testLine[2]))

            # Check non-Nan line, this is the case for any additional header lines
    #         if ~np.isnan(testVals).all()
            if ~np.isnan(testVals[0]):
                EDCS.append(testVals)


    return np.asarray(EDCS), attribs



# Functional form for parsing full set of mat elements and putting in xarray
def EDCSSegsParseX(dumpSegs):
    """
    Extract data from ePS EDCS segments into usable form.

    Parameters
    ----------
    dumpSegs : list
        Set of dumpIdy segments, i.e. dumpSegs, as returned by :py:func:`epsproc.IO.EDCSFileParse()`

    Returns
    -------
    xr.array
        Xarray data array, containing cross sections.
        Dimensions (Eke, theta)

    int
        Number of blank segments found.
        (CURRENTLY not implemented.)


    Example
    -------

    >>> data = EDCSSegsParseX(dumpSegs)

    Notes
    ------

    A rather cut-down version of :py:func:`epsproc.IO.dumpIdySegsParseX()`, no error checking currently implemented.

    """

    dataList = []
    dataArray = []
    ekeList = []
    blankSegs = 0

    # Loop over DumpIdy segments, extract data & reformat
    # If blank, skip parser and append blankSegs.
    for dumpSeg in dumpSegs:
        if len(dumpSeg)>4:
            segBlock, attribs = EDCSSegParse(dumpSeg)
            dataArray.append(segBlock[:,1])         # For brevity, just stack XS data here - will save Xarray sorting later.
            dataList.append([segBlock[:,0], segBlock[:,1], attribs])
            ekeList.append(attribs[0][1])

        else:
            blankSegs += 1
            ekeList.append(np.nan)

    # Dump lists into Xarray - will work provided same theta over all records.
    daOut = xr.DataArray(np.asarray(dataArray), coords={'E':ekeList, 'Theta':segBlock[:,0]}, dims = ['E','Theta'])

    daOut.attrs['dataType'] = 'EDCS'    # Set dataType for use later.

    # Set units - should set from file ideally.
    daOut.attrs['units'] = 'Angs^2'
    daOut.E.attrs['units'] = 'eV'
    daOut.Theta.attrs['units'] = 'deg.'

    return daOut, blankSegs


# ************* CrossSection parsing
# Basically same as EDCS, except read a table of values per symmetry.

# Parse a getCro/CrossSection segment (roughly)
def getCroSegParse(dumpSeg):
    """
    Extract values from GetCro/CrossSection file segments.

    Parameters
    ----------
    dumpSeg : list
        One CrossSection segment, from dumpSegs[], as returned by :py:func:`epsproc.IO.getCroFileParse()`

    Returns
    -------
    np.array
        CrossSections, table of results vs. energy.
    list
        attribs, list [Label, value, units]

    Notes
    -----
    Currently this is a bit messy, and relies on fixed CrossSection output format.
    No error checking as yet.
    Not yet reading all attribs.

    Example
    -------

    >>> XS, attribs = getCroSegParse(dumpSegs[0])

    """
    # Use lists to collect data, and convert format at the end
    # attribs = []
    XS = []

    # attribs.append(['E', np.float(parseLineDigits(dumpSeg[13][2])[0]), 'eV'])
    # attribs.append(dumpSeg[1][2].split())  # Set header line
    attribs = dumpSeg[1][2].split()  # Set header line
    # For each line convert to float - bit ugly, but works
    for testLine in dumpSeg[2:]:
        XS.append(np.genfromtxt(StringIO(testLine[2])))

    return np.asarray(XS), attribs

# Functional form for parsing full set of mat elements and putting in xarray
def getCroSegsParseX(dumpSegs, symSegs, ekeList):
    """
    Extract data from ePS getCro/CrossSecion segments into usable form.

    Parameters
    ----------
    dumpSegs : list
        Set of dumpIdy segments, i.e. dumpSegs, as returned by :py:func:`epsproc.IO.getCroFileParse()`

    Returns
    -------
    xr.array
        Xarray data array, containing cross sections.
        Dimensions (Eke, theta)

    int
        Number of blank segments found.
        (CURRENTLY not implemented.)


    Example
    -------

    >>> data = getCroSegsParseX(dumpSegs)

    Notes
    ------

    A rather cut-down version of :py:func:`epsproc.IO.dumpIdySegsParseX()`, no error checking currently implemented.

    """
    # Old skool debug code.
    # print(ekeList)
    # print(type(ekeList))
    #
    # if (type(ekeList) is not np.ndarray) and (type(ekeList) is not list):
    #     ekeList = [ekeList]   # Wrap as list for single eKE case.
    #     print('Converted to list')

    dataList = []
    dataArray = []
    #ekeList = []
    blankSegs = 0

    # Loop over DumpIdy segments, extract data & reformat
    # If blank, skip parser and append blankSegs.
    for n, dumpSeg in enumerate(dumpSegs):
        if len(dumpSeg)>1:
            segBlock, attribs = getCroSegParse(dumpSeg)

            # Create Xarray - basic
            # daTmp = xr.DataArray(np.asarray(segBlock[:,2:]),
            #         coords={'Ehv':segBlock[:,1], 'XC data':attribs[1:-1:2]}, # 'Type':attribs[2:-1:2]},
            #         dims = ['Ehv', 'XC data']) #, 'Type'])

            # Create Xarray - set MultiIndex first & then rearrange.
            tList = [c[0] for c in attribs[2::2]]  # Set type as char
            typesPD = pd.MultiIndex.from_arrays([attribs[1:-1:2], tList], names = ['XC', 'Type'])
            daTmp = xr.DataArray(np.asarray(segBlock[:,2:]),
                    coords={'Ehv':segBlock[:,1], 'Ctype':typesPD},
                    dims = ['Ehv', 'Ctype']).unstack()

            # Add singleton dim and store.
            # Note this assumes len(symSegs) >= len(dumpSegs)
            # Symmetries are not listed in getCro output.

            sRep = {'ScatSym':'Total', 'ScatContSym':'Cont'}    # Dicitonary look-up for sym names
            try:
                # dataList.append(daTmp.expand_dims({'Sym':symSegs[n]}))
                symList = [symString[2].split('#')[0].split()[1][1:-1] for symString in symSegs[n]]
                symRep = [stringRepMap(symString[2].split()[0],sRep) for symString in symSegs[n]]

            except IndexError as e:
                if e.args[0] != 'list index out of range':
                    raise
                # dataList.append(daTmp.expand_dims({'Sym':'Missing'}))
                symList = ['All', 'All']
                symRep = ['Total', 'Cont']

            # Set MultiIndex for syms
            Syms = pd.MultiIndex.from_arrays([[symList[0]], [symList[1]]], names=symRep)  # Works as expected, rather ugly!
            dataList.append(daTmp.expand_dims({'Sym':Syms}))

        else:
            blankSegs += 1
            #ekeList.append(np.nan)

    # Stack lists by symmetry - this is currently assumed from symSegs
    # if len(dumpSegs) == len(symSegs):
    daOut = xr.combine_nested(dataList, concat_dim = ['Sym'])

    daOut.attrs['dataType'] = 'XSect'    # Set dataType for use later.

    # Set units - should set from file ideally.
    daOut.Ehv.attrs['units'] = 'eV'
    daOut.attrs['units'] = 'Mb'

    # Reset energies to Eke, and shift key dim - might be a simpler/shorter way to do this...?
    # This fails for singleton Ehv/Eke?
    # daOut['EhvOrig'] = daOut['Ehv']
    # daOut['Ehv'] = ekeList
    # daOut = daOut.rename({'Ehv':'Eke', 'EhvOrig':'Ehv'})

    # Try a different approach... assign as new dim then swap.
    if daOut.Ehv.size == 1:
        daOut['Eke'] = ('Ehv', [ekeList])  # Fix for singleton dim - this fails otherwise for scalar or np.ndarray case.
    else:
        daOut['Eke'] = ('Ehv', ekeList)

    daOut = daOut.swap_dims({'Ehv':'Eke'})

    return daOut, blankSegs



# ************* MatEle parsing/sorting - needs a tidy up.

# Linear version of code, for very specific cases.
# Linear tree ip > mu >it
# UGLY... also not working now?  May have passed Xarray dataset in testing by mistake?
def matEleGroupDimX(daIn):
    """
    Group ePS matrix elements by redundant labels (Xarray version).

    Group by ['ip', 'it', 'mu'] terms, all have only a few values.
    Rename 'ip':1,2 as 'Type':'L','V'

    TODO: better ways to do this?
    Via Stack/Unstack? http://xarray.pydata.org/en/stable/api.html#id16
    See also tests in funcTests_210819.py for more versions/tests.

    Parameters
    ----------
    data : Xarray
        Data array with matrix elements to be split and recombined by dims.

    Returns
    -------
    data : Xarray
        Data array with reordered matrix elements (dimensions).

    """

    daRedList = []
    daRed = daIn

    # Split on mu - values from set {-1,0,1} depending on symmetry
    daRedList = []
    uVals = np.unique(daRed.mu)
    for n, val in enumerate(uVals):
        tmp = matEleSelector(daRed, inds = {'mu':val})
        tmp = tmp.expand_dims({'mu':[val]})
        daRedList.append(tmp)

    # Restack
    daRed = xr.combine_nested(daRedList, concat_dim = 'mu')

    # Split on it
    daRedList = []
    uVals = np.unique(daRed.it)
    for n, val in enumerate(uVals):
        tmp = matEleSelector(daRed, inds = {'it':val})
        tmp = tmp.expand_dims({'it':[val]})
        daRedList.append(tmp)

    # Restack
    daRed = xr.combine_nested(daRedList, concat_dim = 'it')

    # Split on 'ip' - will always be (1,2), and split matEle into two
    ipLabel = ['L','V']
    daRedList = []
    for n, val in enumerate(range(1,3)):
        tmp = matEleSelector(daRed, inds = {'ip':val})
        tmp = tmp.expand_dims({'Type':[ipLabel[n]]})
        daRedList.append(tmp)

    # Restack
    daRed = xr.combine_nested(daRedList, concat_dim = 'Type')

    return daRed

#   Subselections using matEleSelector
#   THIS IS UGLY, but seems to work consistently - get da of correct dims out (with multilevel coords in).
#   Could also try dataset to array for split and recombine...?
#   http://xarray.pydata.org/en/v0.12.3/reshaping.html
def matEleGroupDimXnested(da):
    """
    Group ePS matrix elements by redundant labels (Xarray version).

    Group by ['ip', 'it', 'mu'] terms, all have only a few values.

    TODO: better ways to do this?
    See also tests in funcTests_210819.py for more versions/tests.

    Parameters
    ----------
    data : Xarray
        Data array with matrix elements to be split and recombined by dims.

    """

    indList = ['ip','it','mu']

    daRedList = []
    for x in np.unique(da[indList[0]]):
        daRedList0 = []
        for y in np.unique(da[indList[1]]):
            daRedList1 = []
            for z in np.unique(da[indList[2]]):
                red = matEleSelector(da, inds = {indList[0]:x, indList[1]:y, indList[2]:z})

                red = red.expand_dims({indList[0]:[x], indList[1]:[y], indList[2]:[z]})

                daRedList1.append(red)

            daOut1 = xr.combine_nested(daRedList1, concat_dim = indList[2])
            daRedList0.append(daOut1)

        daOut2 = xr.combine_nested(daRedList0, concat_dim = indList[1])
        daRedList.append(daOut2)

    daOut =  xr.combine_nested(daRedList, concat_dim = indList[0])

    return daOut

# UGH THIS IS SO UGLY, please make it better.
# Should be a neat recursive tree method here, probably also something canned!
# Or with native Xarray functionality, but in testing couldn't get this to work properly.
def matEleGroupDim(data, dimGroups = [3, 4, 2]):
    """
    Group ePS matrix elements by redundant labels.

    Default is to group by ['ip', 'it', 'mu'] terms, all have only a few values.

    TODO: better ways to do this?  Shoud be possible at Xarray level.

    Parameters
    ----------
    data : list
        Sections from dumpIdy segment, as created in dumpIdySegsParseX()
        Ordering is [labels, matElements, attribs].

    """

    # Basic version assuming dims

#    # Split on ip
    dataIP = []
    dInd = 3
    dataIP = dataGroupSel(data, dInd)

    # Split on IT
    dInd = 4
    dataIT = []
    for dataSub in dataIP:
        temp = dataSub
        dataIT.extend(dataGroupSel(dataSub, dInd))


    # Split on mu
    dInd = 2
    dataMU = []
    for dataSub in dataIT:
        dataMU.extend(dataGroupSel(dataSub, dInd))


    # Put into Xarray (code from dumpIdySegsParseX)
    # Label singleton dims directly, then stack.
    dataArrays = []
    attribs = data[2]   # Shared attribs
    for dataSub in dataMU:
        LM = pd.MultiIndex.from_arrays(dataSub[0][0:2,:].astype('int8'), names = attribs[-1][1][0:2])
        LM = LM.swaplevel(0, 1)  # Switch l,m indexes
        mu = [dataSub[0][2,0].astype('int8')]    # Already set to single values above.
        ip = [dataSub[0][3,0].astype('int8')]
        it = [dataSub[0][4,0].astype('int8')]
        Syms = pd.MultiIndex.from_tuples([(attribs[4][1],attribs[5][1],attribs[6][1])],names=[attribs[4][0],attribs[5][0],attribs[6][0]])

        #dataArrays.append(xr.DataArray(data[1], coords={'ES': Esyms, 'LM':QNs}, dims = ['ES','LM']))
        # AH - issue is number of labels - can't lable singleton dim it seems, but can expand
        #TODO: consider setting E as a separate dim, will be singleton for each set of syms. Might make more sense for later manipulations (sum over sym or E).
        # tmp = xr.DataArray(np.asarray(data[1]), coords={'LM':QNs}, dims = ['LM'])
        # tmp = tmp.expand_dims({'Sym':Syms, 'Eke':[attribs[0][1]]})  # This is OK, but still ties Eke and Sym coords (same number of elements)
        # # tmp = tmp.expand_dims({'Sym':Syms})

        # tmp =  xr.DataArray(np.asarray(data[1]), coords={'LM':LM, 'mu':mu, 'ip':ip, 'it':it}, dims = ['LM','mu','it','ip'])
        tmp =  xr.DataArray(np.asarray(dataSub[1]), coords={'LM':LM}, dims = ['LM'])
        tmp = tmp.expand_dims({'Sym':Syms, 'Eke':[attribs[0][1]]})
        tmp = tmp.expand_dims({'mu':mu, 'ip':ip, 'it':it})

         # Assign any other attributes - note that some attributes may be dropped when combining arrays below
        for a in attribs:
            tmp.attrs[a[0]] = a[1] # Currently set without units, multiple values here give combine issues below.

        dataArrays.append(tmp)

    # Recombine along it
    da = xr.combine_nested(dataArrays, concat_dim = ['it'])

    return da




# Function for grabbing files or scanning dir for files.
# Note raw string for docstring as one method of keeping raw string in example.
def getFiles(fileIn = None, fileBase = None, fType = '.out', verbose = True):
    r"""
    Read ePS file(s) and return results as Xarray data structures.
    File endings specified by fType, default .out.

    Parameters
    ----------
    fileIn : str, list of strs, optional.
        File(s) to read (file in working dir, or full path).
        Defaults to current working dir if only a file name is supplied.
        For consistent results, pass raw strings, e.g.
        ``fileIn = r"C:\share\code\ePSproc\python_dev\no2_demo_ePS.out"``

    fileBase : str, optional.
        Dir to scan for files.
        Currently only accepts a single dir.
        Defaults to current working dir if no other parameters are passed.

    fType : str, optional
        File ending for ePS output files, default '.out'

    verbose : bool, optional
        Print output details, default True.


    Returns
    -------
    list
        List of Xarray data arrays, containing matrix elements etc. from each file scanned.


    """


    currDir = os.getcwd()

    if fileBase is None:
        fileBase = currDir

    if fileIn is not None:
        # Wrap in list if only single file passed
        if type(fileIn) is not list:
            fileIn = [fileIn]

        fList = []
        for file in fileIn:
            # Check file & path are valid
            fTest = os.path.split(file)
            if not fTest[0]:
                fList.append(os.path.join(currDir, file))
            else:
                fList.append(file)

        # Display message
        if verbose:
            print('\n*** Scanning file(s)')
            print(fList)

    else:
        # Filenames only
        # fList = [f for f in os.listdir(fileBase) if f.endswith(fType)]
        # With full path
        fList = [os.path.join(fileBase, f) for f in os.listdir(fileBase) if f.endswith(fType)]

        # Display message
        if verbose:
            print('\n*** Scanning dir')
            print(fileBase)
            print('Found {0} {1} file(s)\n'.format(len(fList), fType))

    return fList

#****** Master function to read a file, or dir, of ePS outputs.

# Some of the logic and methods could do with a revisit/tidy-up here...
#TODO: Add error checking on paths using os.path.isdir/isfile etc.
#TODO: Check/fix paths if incorrectly passed, e.g. https://stackoverflow.com/a/21605790
# ADDED: type switch for matEle or EDCS, probably more to come. Should rename function!

def readMatEle(fileIn = None, fileBase = None, fType = '.out', recordType = 'DumpIdy', verbose = 1, stackE = True):
    r"""
    Read ePS file(s) and return results as Xarray data structures.
    File endings specified by fType, default *.out.

    Parameters
    ----------
    fileIn : str, list of strs, optional.
        File(s) to read (file in working dir, or full path).
        Defaults to current working dir if only a file name is supplied.
        For consistent results, pass raw strings, e.g.
        ``fileIn = r"C:\share\code\ePSproc\python_dev\no2_demo_ePS.out"``

    fileBase : str, optional.
        Dir to scan for files.
        Currently only accepts a single dir.
        Defaults to current working dir if no other parameters are passed.

    fType : str, optional
        File ending for ePS output files, default '.out'

    recordType : str, optional, default 'DumpIdy'
        Type of record to scan for, currently set for 'DumpIdy', 'EDCS' or 'CrossSection'.
        For a full list of descriptions, types and sources, run:
        >>> epsproc.util.dataTypesList()

    verbose : int, optional, default = 1
        Level of verbosity in output.
        - 0 no printed output
        - 1 print summary info only
        - 2 print detailed info

    stackE : bool, optional, default = True
        Identify and stack multi-part jobs to single array (by E) if True.

    Returns
    -------
    list
        List of Xarray data arrays, containing matrix elements etc. from each file scanned.

    To do
    -----
    - Change to pathlib paths.
    - Implement outputType options...?

    13/10/20    Adapted to use grouped lists for multi-file jobs, should be back-compatible if stackE = False set.

    Examples
    --------

    >>> dataSet = readMatEle()  # Scan current dir

    >>> fileIn = r'C:\share\code\ePSproc\python_dev\no2_demo_ePS.out'
    >>> dataSet = readMatEle(fileIn)  # Scan single file

    >>> dataSet = readMatEle(fileBase = r'C:\share\code\ePSproc\python_dev') # Scan dir

    .. note::
        * Files are scanned for matrix element output flagged by "DumpIdy" headers.
        * Each segment found is parsed for attributes and data (set of matrix elements).
        * Matrix elements and attributes are combined and output as an Xarray array.

    """

    # Set files to read, either:
    #   - dir (one only) to scan using os
    #   - file(s) and add to list with full path
    #   - default case is scan current working dir.

    if verbose:
        print('*** ePSproc readMatEle(): scanning files for ' + recordType + ' segments.')

    # Call function to check files or scan dir.
    fList = getFiles(fileIn = fileIn, fileBase = fileBase, fType = fType, verbose = verbose)

    # Check if job is multipart, and stack files if so, based on filename prefix
    if stackE:
        fListSorted, fListGrouped, prefixStr = fileListSort(fList, verbose = verbose)

    else:
        fListGrouped = [fList]  # Default case, just use fList here.
        prefixStr = ''

    if verbose > 2:
        print(fListGrouped)


    # Loop over fList and scan ePS files
    # 13/10/20 Adapted to use grouped subsets.
    dataSet = []
    dataStack = []
    for fList in fListGrouped:
        # Force to list for singleton case.
        if not isinstance(fList, list):
            fList = [fList]

        # Check if files should be Eke stacked
        if stackE and len(fList) > 1:
            stackFlag = True
        else:
            stackFlag = False

        for file in fList:
            if verbose:
                print('\n*** Reading ePS output file: ', file)

            # Scan the file and parse segments
            #lines, dumpSegs = dumpIdyFileParse(os.path.join(fileBase, file))

            if recordType == 'DumpIdy':
                ekeList = scatEngFileParse(file, verbose = verbose)
                symSegs = symFileParse(file, verbose = verbose)

                if verbose > 1:
                    print('Scanning CrossSection segments.')
                    print('Expecting {0} DumpIdy segments.'.format(ekeList.size * len(symSegs)))

                lines, dumpSegs = dumpIdyFileParse(file, verbose = verbose)
                data, blankSegs = dumpIdySegsParseX(dumpSegs, ekeList, symSegs, verbose = verbose)

            if recordType == 'EDCS':
                # print('Expecting {0} EDCS segments.'.format(ekeList.size))
                if verbose > 1:
                    print('Scanning EDCS segments.')

                lines, dumpSegs = EDCSFileParse(file, verbose = verbose)
                data, blankSegs = EDCSSegsParseX(dumpSegs) # , ekeList, symSegs)

            if recordType == 'CrossSection':
                ekeList = scatEngFileParse(file, verbose = verbose)
                symSegs = symFileParse(file, verbose = verbose)

                if verbose > 1:
                    print('Scanning CrossSection segments.')
                    print('Expecting {0} CrossSection segments.'.format(len(symSegs)+1))  # Assuming 1 segment per symmetry, plus symm-summed case.

                lines, dumpSegs = getCroFileParse(file, verbose = verbose)
                data, blankSegs = getCroSegsParseX(dumpSegs, symSegs, ekeList)



            # Add some additional properties to the output
            fName = os.path.split(file)
            data.name = fName[1]
            data.attrs['file'] = fName[1]
            data.attrs['fileBase'] = fName[0]

            if verbose:
                print('Processed {0} sets of {1} file segments, ({2} blank)'.format(len(dumpSegs),recordType,blankSegs))

            # For unstaked case, append each file to dataSet
            if stackFlag:
                dataStack.append(data)
            else:
                dataSet.append(data)

        # Stack by Eke
        if stackFlag:
            dataSet.append(xr.combine_nested(dataStack, concat_dim = ['Eke']).sortby('Eke'))

            # Propagate attribs for stacked case
            dataSet[-1].attrs = dataStack[0].attrs
            dataSet[-1].attrs['fileList'] = fList

            if verbose:
                print(f"\n*** Stacked {len(fList)} files, prefix {prefixStr}, by Eke ({dataSet[-1].Eke.size} points).")

        else:
            dataSet[-1].attrs['fileList'] = dataSet[-1].attrs['file']  # Set also for single file case

        # if outputType == 'list':
        # Put in a list for now, might want to use Xarray dataset here, and/or combine results from multiple files.

        # if outputType == 'ds':
            # Output to Xarray Dataset, including stacking by E for multi-file jobs.



    return dataSet


# **************** Functions for wavefunction (3D data) files
# Based on previous Matlab code, readOrb3D.m
# Helped by Scipy Cookbook LAS reader example: https://scipy-cookbook.readthedocs.io/items/LASReader.html
#
# See:
#   - ePSproc_dev_IOfun_260519.py
#   - ePSproc_dev_3Dvis_290519.property
# for development notes (Finn E:\ePS_paraview_proc\scripts).
#
# TODO: rename & integrate with functions above.

# Read header lines into a list & convert to int
def readOrbHeader(f):
    headerLines=[]
    for i in range(5):
        headerLines.append(f.readline())
        # headerLines.append(f.readline().split()) # Split at whitespace

        if i>0:
            headerLines[i] = int(headerLines[i])  # Convert to int

    return headerLines

# Read a specified number of floats from file, return as numpy array type
def readOrbElements(f,n):
    data = []

    while len(data) < n:
        data.extend([float(s) for s in f.readline().split()])

    return np.array(data)

# Read coords from file, based on params in headerLines
def readOrbCoords(f, headerLines):
    coords = []
    for n in range(np.abs(headerLines[1])):
        coords.append(readOrbElements(f,headerLines[n+2]))

    return coords

# Read data from file, based on params in headerLines
def readOrbData(f, headerLines):

    # Determine size of grid
    nData = 1
    nGrid = []
    # [j=j*headerLines[k] for k in range(np.abs(headerLines[1]))]  # FAIL
    # Set number of elements
    for n in range(np.abs(headerLines[1])):
        nData = nData*headerLines[n+2]
        nGrid.append(headerLines[n+2])


    # Read elements
    dataRaw = readOrbElements(f,nData)

    # Resort to 3D array - adapted from Matlab code, so might be a +/-1 index offset mix-up here...
    C = 0
    data = np.zeros(nGrid)
    for z in range(nGrid[2]):
        for y in range(nGrid[1]):
           # data[0:(nGrid[0]-1),y,z] = dataRaw[C:(C+nGrid[0]-1)]  # Incorrect offset?
           data[:,y,z] = dataRaw[C:(C+nGrid[0])]  # Should be correct... implicit (-1) in dataRaw range setting
           C = C+nGrid[0]


    return data

# *************** Master function for reading a set of 3D data files from ePS
def readOrb3D(fileIn = None, fileBase = None, fType = '_Orb.dat', verbose = True):
    """
    Read ePS 3D data file(s) and return results.
    File endings specified by fType, default *_Orb.dat.

    Parameters
    ----------
    fileIn : str, list of strs, optional.
        File(s) to read (file in working dir, or full path).
        Defaults to current working dir if only a file name is supplied.
        For consistent results, pass raw strings, e.g.
        fileIn = r"C:\share\code\ePSproc\python_dev\no2_demo_ePS.out"

    fileBase : str, optional.
        Dir to scan for files.
        Currently only accepts a single dir.
        Defaults to current working dir if no other parameters are passed.

    fType : str, optional
        File ending for ePS output files, default '_Orb.dat'

    verbose : bool, optional
        Print output details, default True.


    Returns
    -------
    list
        List of data arrays, containing matrix elements etc. from each file scanned.

    # TODO: Change output to Xarray?

    Examples
    --------

    >>> dataSet = readOrb3D()  # Scan current dir

    >>> fileIn = r'C:\share\code\ePSproc\python_dev\DABCOSA2PPCA2PP_10.5eV_Orb.dat'
    >>> dataSet = readOrb3D(fileIn)  # Scan single file

    >>> dataSet = readOrb3D(fileBase = r'C:\share\code\ePSproc\python_dev') # Scan dir


    """

    # Populate file list
    fList = getFiles(fileIn = fileIn, fileBase = fileBase, fType = fType, verbose = verbose)

    dataSet = []
    for fileName in fList:
        with open(fileName,'r') as f:
            # Check eof
            f.seek(0,2)
            fEnd = f.tell()
            f.seek(0)

            # Read file segments until eof
            # TODO: add eof checking here, can have 3 or 6 segments depending on symmetry.
            data = []
            # for seg in range(3):
            while f.tell() < fEnd - 100:    # Scan until eof minus arb small offset.
                headerLines = readOrbHeader(f)
                coords = readOrbCoords(f,headerLines)
                data.append(readOrbData(f, headerLines))

        dataSet.append([fileName, headerLines, coords, data])

    return dataSet


def writeOrb3Dvtk(dataSet):
    """
    Write ePS 3D data file(s) to vtk format.
    This can be opened in, e.g., Paraview.

    Parameters
    ----------
    dataSet : list
        List of data arrays, containing matrix elements etc. from each file scanned.
        Assumes format as output by readOrb3D(), [fileName, headerLines, coords, data]

    # TODO: Change to Xarray?

    Returns
    -------
    list
        List of output files.


    .. note::
        Uses Paulo Herrera's eVTK, see:

        * https://pyscience.wordpress.com/2014/09/06/numpy-to-vtk-converting-your-numpy-arrays-to-vtk-arrays-and-files/
        * https://bitbucket.org/pauloh/pyevtk/src/default/

        `pip install pyevtk` to install.


    """
    fOut = []
    for file in dataSet:
        # Set grid, convert to Cart if necessary, assuming that grid won't be larger than 10 Angs
        # if (len(file[2][0]) > 50):
        #     # Convert to Cart grid for plotting
        #     # TODO: Investigate use of sph grid here - should be cleaner.
        #     # TODO: Investigate recreating mesh in Paraview, rather than saving to file.
        #     [T,R,P] = np.meshgrid(file[2][1], file[2][0], file[2][2])
        #     T = (T*np.pi/180) #-np.pi/2
        #     P = P*np.pi/180
        #     x = R*np.sin(P)*np.cos(T)
        #     z = R*np.cos(P)
        #     y = R*np.sin(P)*np.sin(T)
        # else:
        #     [x,y,z] = np.meshgrid(file[2][1], file[2][0], file[2][2])

        # Now set as separate function
        [x,y,z] = orb3DCoordConv(file)

        # Save single dataset
        # Function info: https://bitbucket.org/pauloh/pyevtk/src/default/src/hl.py
        # gridToVTK("./ePStest3", x, y, z, pointData = {'data': data[0]})

        # Save full dataset
        # TODO: Check number of segments and save multiple
        # segs = round(len(file[3])/3)
        fOut.append(gridToVTK(file[0][:-4], x, y, z, pointData = {'Re': file[3][0], 'Im': file[3][1], 'Abs': file[3][2]}))

    print("{0} files written to vtk format.".format(len(fOut)))

    return fOut


#**************** Wrappers for Xarray load/save netCDF

# File write wrapper.
def writeXarray(dataIn, fileName = None, filePath = None, engine = 'h5netcdf', forceComplex = False):
    """
    Write file to netCDF format via Xarray .to_netcdf() method.

    Parameters
    -----------
    dataIn : Xarray
        Data array to write to disk.

    fileName : str, optional, default = None
        Filename to use.
        If set to None (default) the file will be written with a datastamp.

    filePath : str, optional, default = None
        Full path to file.
        If set to None (default) the file will be written in the current working directory (as returned by `os.getcwd()`).

    engine : str, optional, default = 'h5netcdf'
        netCDF engine for Xarray to_netcdf method. Some libraries may not support multidim data formats.
        See https://docs.xarray.dev/en/latest/user-guide/io.html

    forceComplex : bool, optional, default = False
        For h5netcdf engine only, set `invalid_netcdf` option = forceComplex.
        If True, complex data will be written directly to file.
        Note this also needs to be read back with the same engine & settings.
        For more details see https://github.com/h5netcdf/h5netcdf#invalid-netcdf-files

    Returns
    -------
    str
        Indicates save type and file path.

    Notes
    -----
    The default option for Xarray is to use Scipy netCDF writer, which does not support complex datatypes. In this case, the data array is written as a dataset with a real and imag component.

    This routine assumes a DataArray as input, although the Xarray file writer pushes this to a DataSet.

    TODO: implement try/except to handle various cases here, and test other netCDF writers (see http://xarray.pydata.org/en/stable/io.html#netcdf).

    Multi-level indexing is also not supported, and must be serialized first. Ugh.

    02/06/22: added improved complex number handling & attibutes sanitizer (lossy).

    """

    if fileName is None:
        timeString = dt.now()
        fileName = 'ep_' + timeString.strftime('%Y-%m-%d_%H-%M-%S')

    if filePath is None:
        filePath = os.getcwd()

    # Serialize MultiIndex  - testing here for BLM case.
    # if 'BLM' in dataIn.dims:
    #     dataIn = dataIn.reset_index(['Euler','BLM'])

    # Serialize general - use unstact() to flatten all dims
    dataIn = dataIn.unstack()

# Try/except not yet working, multiple error types to handle here...
    # try:
    #     dataIn.to_netcdf(fileName)
    #     saveMsg = 'Written to netCDF4.'
    #     print(saveMsg)
    #     return saveMsg
    #
    # except ValueError as e:
    #     if e.msg != "NetCDF 3 does not support type complex128":
    #         raise
    #     else:
    #         xr.Dataset({'Re':dataIn.real, 'Im':dataIn.imag}).to_netcdf(fileName)
    #         saveMsg = 'Written to netCDF3 (re/im format).'
    #         print(saveMsg)
    #         return saveMsg
    #
    # return 'File not written.'


    if (engine != 'h5netcdf') or (not forceComplex):
        # Safe version with re/im split save type only.
        # Works for scipy and h5netcdf OK, latter will save complex type too, but is strictly not valid.
        dataOut = xr.Dataset({'Re':dataIn.real, 'Im':dataIn.imag})
        # dataOut.attrs = dataIn.attrs   # This will push dataarray attrs to dataset attrs, otherwise they're nested
                                        # May not always want this?

        # Allow for SF & XS coords which may also be complex
        # if 'XS' in dataOut.coords:
        #     dataOut['XSr'] = dataOut.XS.real
        #     dataOut['XSi'] = dataOut.XS.imag
        #     dataOut = dataOut.drop('XS')
        #
        # if 'SF' in dataOut.coords:
        #     dataOut['SFr'] = dataOut.SF.real
        #     dataOut['SFi'] = dataOut.SF.imag
        #     dataOut = dataOut.drop('SF')

        # Allow for arb complex coords.
        # May also want to add attr checker here? Or set in 'sanitizeAttrsNetCDF'
        for item in dataOut.coords.keys():
            if dataOut.coords[item].dtype == 'complex128':
                dataOut.coords[item + 'r'], dataOut.coords[item + 'i'] = splitComplex(dataOut.coords[item])
                dataOut = dataOut.drop(item)

    else:
        # dataOut = dataIn.to_dataset()   # Set to dataset explicitly prior to save - may also need/want to set name here if missing.
        dataOut = dataIn   # Without additional conversion.

    # For netCDF3 can't have multidim attrs, quick fix here for removing them (BLM case)
    if engine == 'scipy':
        if 'sumDims' in dataOut.attrs:
            dataOut.attrs['sumDims'] = [] # test.attrs['selDims'][0]
        if 'selDims' in dataOut.attrs:
            dataOut.attrs['selDims'] = []

    try:
        if engine != 'h5netcdf':
            dataOut.to_netcdf(os.path.join(filePath, fileName + '.nc'), engine=engine)
        else:
            dataOut.to_netcdf(os.path.join(filePath, fileName + '.nc'), engine=engine, invalid_netcdf=forceComplex)

        saveMsg = [f'Written to {engine} format']

    except Exception as e:

        print(f'writeXarray caught exception: {e}')
        print(f'Retrying file write with sanitized attrs.')

        # THIS IS WORKING IN TESTING, but not here?
        # Does work for invalid_netcdf = True case, for h5netcdf backend at least.
        # Seems to be an issue with sanitizeAttrsNetCDF() and/or backend to fix?
        # AH - issue is DataSet vs. DataArray. In former case may need to fix attrs per data variable too.
        # TODO: make use of sanitizeAttrsNetCDF log return. Write to sidecar file?
        dataOut, attrs, log = sanitizeAttrsNetCDF(dataOut)
        if isinstance(dataOut, xr.core.dataset.Dataset):
            for item in dataOut.data_vars:
                dataOut[item], attrs, log = sanitizeAttrsNetCDF(dataOut[item])

        if engine != 'h5netcdf':
            dataOut.to_netcdf(os.path.join(filePath, fileName + '.nc'), engine=engine)
        else:
            dataOut.to_netcdf(os.path.join(filePath, fileName + '.nc'), engine=engine, invalid_netcdf=forceComplex)

        saveMsg = [f'Written to {engine} format, with sanitized attribs (may be lossy)']

    saveMsg.append(os.path.join(filePath, fileName + '.nc'))
    print(saveMsg)

    return saveMsg


# Split complex to R + I
def splitComplex(data):
    """Split complex data into R+I floats."""

    dataR = np.real(data)
    dataI = np.imag(data)

    return dataR, dataI

# Comibine R + I to complex
def combineComplex(dataR, dataI):
    """Combine R+I floats into complex form."""

    data = dataR + 1j*dataI

    return data


# Sanitize attributes & dicts for Xarray NetCDF IO
def sanitizeAttrsNetCDF(data, dictHandling = 'wrap'):
    """
    Sanitize Xarray DataArray attributes for file IO.

    Note this may be lossy:

    - Empty data > string.
    - Dictionaries removed, wrapped to string, or left alone (nested dicts not supported in attrs for most (all?) file writers).
      Set dictHandling = 'del', 'wrap' or anything else to leave as is.
    - Remove all items not of types [str, np.ndarray, int, float, list, tuple]


    Todo:

    - try conversion to string for all attrs?
    - try dict conversions & JSON side-car file IO to avoid lossy saves.

    """

    dataOut = data.copy()

    # Remove None and other empty types, ugh - now integrated below
    # xrTest.attrs = {k:(v if v else str(v)) for k,v in xrTest.attrs.items()}
    log = {}
    for k,v in dataOut.attrs.items():
        if not v:
            dataOut.attrs[k] = str(v)
            log[k] = 'str'

        if isinstance(dataOut.attrs[k], dict):
    #         xrTest.attrs[k] = [[k2,v2] for k2,v2 in xrTest.attrs[k].items()]  # Nest dict items also not supported, dump to nested lists? Seems to be acceptable. Ugh.
                                                                                # Still causing issues in some cases?
            if dictHandling == 'del':
                dataOut.attrs[k] = 'Removed dict'
                log[k] = 'Removed dict'

            elif dictHandling == 'wrap':
                dataOut.attrs[k] = str(v)
                log[k] = 'Wrapped dict to string'

            else:
                pass

        if type(dataOut.attrs[k]) not in [str, np.ndarray, int, float, list, tuple]:
            typeIn = type(dataOut.attrs[k])
            dataOut.attrs[k] = 'NA'
            log[k] = f'Removed item type {typeIn}'

    # TO TRY - full str conversion, e.g. from https://stackoverflow.com/a/42676094 (for JSON example case)
    # save: convert each tuple key to a string before saving as json object
    # s = json.dumps({str(k): str(v) for k, v in eulerDict.items()})
    #
    # THEN RECON with ast:
    # # load in two stages:
    # # (i) load json object
    # obj = json.loads(s)
    #
    # # (ii) convert loaded keys from string back to tuple
    # from ast import literal_eval
    # # d = {literal_eval(k): literal_eval(v) for k, v in obj.items()}  # FAILS: ValueError: malformed node or string: <ast.Name object at 0x7f4464e67550>
    # d = {k: (literal_eval(v) if v != 'Euler' else v) for k, v in obj.items()}  # ok - WORKS FOR ALL CASES EXCEPT NON-EXECUTABLE STRS
    #
    # This should also work here, but maybe add type checking too?


    return dataOut, data.attrs, log


# File read wrapper.
def readXarray(fileName, filePath = None, engine = 'h5netcdf', forceComplex = False, forceArray = True):
    """
    Read file from netCDF format via Xarray method.

    Parameters
    -----------
    fileName : str
        File to read.

    filePath : str, optional, default = None
        Full path to file.
        If set to None (default) the file will be written in the current working directory (as returned by `os.getcwd()`).


    engine : str, optional, default = 'h5netcdf'
        netCDF engine for Xarray to_netcdf method. Some libraries may not support multidim data formats.
        See https://docs.xarray.dev/en/latest/user-guide/io.html

    forceComplex : bool, optional, default = False
        For h5netcdf engine only, set `invalid_netcdf` option = forceComplex.
        If True, complex data will be written directly to file.
        Note this also needs to be read back with the same engine & settings.
        For more details see https://github.com/h5netcdf/h5netcdf#invalid-netcdf-files

    forceArray : bool, optional, default = False
        Force file reader to use xr.open_dataarray if True.
        Otherwise use xr.open_dataset.
        The latter case may need additional post-processing, but works with cases with split Re & Im dataarrays.
        If forceComplex = False this setting is ignored.


    Returns
    -------
    Xarray
        Data from file.  May be in serialized format.

    Notes
    -----
    The default option for Xarray is to use Scipy netCDF writer, which does not support complex datatypes. In this case, the data array is written as a dataset with a real and imag component.

    Multi-level indexing is also not supported, and must be serialized first. Ugh.

    TODO: generalize multi-level indexing here.


    07/06/22: improved dim restacking with :py:func:`epsproc.util.misc.restack` routine.
    02/06/22: improved engine & complex number handling (as per writeXarray)

    """

    # TODO - file and path checks
    # See writeXarray() above, and PEMtk.fit._io() functions - should unify method here!
    # See also ep.IO.getFiles?
    # if dataPath is None:
    #     # dataPath = os.getcwd()  # OR
    #     dataPath = Path().absolute()
    # if not Path(fileIn).exists():
    #     fileIn = Path(dataPath,fileIn)  # Assume full path missing if file doesn't exist?

    # Set reader - can try and force to array too.
    # If forceComplex = False, need to use xr.open_dataset for Re+Im dataset format.
    if forceArray and forceComplex:
        freader = xr.open_dataarray
    else:
        freader = xr.open_dataset

    # Read file
    if engine != 'h5netcdf':
        dataIn = freader(fileName, engine = engine)
    else:
        dataIn = freader(fileName, engine = engine, invalid_netcdf = forceComplex)

    if (engine != 'h5netcdf') or (not forceComplex):
        # Reconstruct complex variables, NOTE this drops attrs... there's likely a better way to do this!
        # UPDATE 07/06/22: additional attrs handling below. Note in this case dataOut is a DataArray here.
        dataOut = dataIn.Re + dataIn.Im*1j
        # dataOut.attrs = dataIn.attrs

        # Rest SF & XS coords which may also be complex
        # Note: need to check vs. dataIn here, since dataOut already has dropped vars
        # if 'XSr' in dataIn.data_vars:
        #     dataOut['XS'] = dataIn.XSr + dataIn.XSi*1j
        # #     dataOut = dataOut.drop('XSr').drop('XSi')
        #
        # if 'SFr' in dataIn.data_vars:
        #     dataOut['SF'] = dataIn.SFr + dataIn.SFi
        # #     dataOut = dataOut.drop('SFr').drop('SFi')

        # General version
        for item in dataOut.coords.keys():
            # Check for r+i pairs - note labelling assumed to match writeXarray conventions here.
            if item.endswith('r'):
                itemi = item[:-1] + 'i'

                # If imag partner found, restack and remove split components.
                if itemi in dataOut.coords.keys():
                    dataOut.coords[item[:-1]] = combineComplex(dataOut.coords[item], dataOut.coords[itemi])
                    dataOut = dataOut.drop([item,itemi])

    else:
        dataOut = dataIn

    # For dataset case, try some generic handling. May need more sophisticated methods here, maybe just assume DataArray and convert?
    if (not dataOut.attrs) and isinstance(dataIn, xr.core.dataset.Dataset):
        dataOut.attrs = dataIn[list(dataIn.data_vars)[0]].attrs

    # Recreate MultiIndex from serialized version  - testing here for BLM case.
    # if 'BLM' in dataIn.dims:
    #     dataIn = dataIn.set_index({'BLM':['l','m'],'Euler':['P','T','C']})

    # Recreate MultiIndex from serialized version according to array type.
    # 01/06/22: added try/except for lazy dim handling.
    try:
        # if dataIn.dataType == 'BLM':
        #     dataOut = dataOut.stack(BLMdimList(sType = 'sDict'))
        # elif dataIn.dataType == 'matE':
        #     dataOut = dataOut.stack(matEdimList(sType = 'sDict'))

        dataOut, dims = restack(dataOut)  # General restacking routine, may want to pass args here for more flexibility.

    except:
        print(f"Failed to restack input dataset for dataType {dataIn.dataType}, dims may be missing. Check ep.dataTypesList['{dataIn.dataType}'] for details.")


    return dataOut
