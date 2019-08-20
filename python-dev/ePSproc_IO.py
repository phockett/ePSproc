# -*- coding: utf-8 -*-
"""
ePSproc IO functions.

Main function readMatEle(fileIn = None, fileBase = None, fType = '.out'):

    Read ePS file(s) and return results as Xarray data structures containing matrix elements.
    File endings specified by fType, default *.out.

Also includes required ancillary functions for file IO tasks, and data parsing.

19/08/19        Add functions for reading wavefunction files (3D data)
07/08/19        Naming convention tweaks, and some changes to comments, following basic tests with Sphinx.
05/08/19    v1  Initial python version.
                Working, but little error checking as yet. Needs some tidying.

TODO
----
- Add IO for other file segments (only DumpIdy supported so far).
- Better logic & flexibility for file scanning.
- Restructure as class for brevity...?

"""

# Imports
import os
import re
import numpy as np
import pandas as pd
from io import StringIO
import xarray as xr

try:
    from pyevtk.hl import gridToVTK
except ImportError as e:
    if e.msg != "No module named 'pyevtk'":
        raise
    print('* pyevtk not found, VTK export not available. ')


# ***** Ancillary functions

# File parsing function - scan file for keywords & read segments
#   Following above idiomatic solution, with full IO
#   https://stackoverflow.com/questions/3961265/get-line-number-of-certain-phrase-in-file-python
def fileParse(fileName, startPhrase = None, endPhrase = None, comment = None):
    """
    Parse a file, return segment(s) from startPhrase:endPhase, excluding comments.

    Parameters
    ----------
    fileName : str
        File to read (file in working dir, or full path)
    startPhrase : str, optional
        Phrase denoting start of section to read. Default = None
    endPhase : str, optional
        Phrase denoting end of section to read. Default = None
    comment : str, optional
        Phrase denoting comment lines, which are skipped. Default = None

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

    # Open file & scan each line.
    with open(fileName,'r') as f:
        for (i, line) in enumerate(f):  # Note enumerate() here gives lines with numbers, e.g. fullFile=enumerate(f) will read in file with numbers

            # If line matches startPhrase, print line & append to list.
            if startPhrase in line:
               print('Found "', startPhrase, '" at line: ', i)
               lineStart.append(i)

               readFlag = True

            # Read lines into segment[] until endPhrase found
            if readFlag:
                # Check for end of segment (start of next Command sequence)
                if endPhrase and line.startswith(endPhrase):
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


    return ([lineStart, lineStop], segments[:-1])


# Simple wrapper for general fileParse function, ePS dumpIdy segments
def dumpIdyFileParse(fileName):
    """
    Parse an ePS file for dumpIdy segments.

    Parameters
    ----------
    fileName : str
        File to read (file in working dir, or full path)

    Returns
    -------
    list
        [lineStart, lineStop], ints for line #s found from start and end phrases.
    list
        dumpSegs, list of lines read from file.

    Lists contain entries for each dumpIdy segment found in the file.

    """

    startPhrase = "DumpIdy - dump"
    endPhrase = "+ Command"

    (lines, dumpSegs) = fileParse(fileName, startPhrase, endPhrase) # , '>')

    return lines, dumpSegs


# Parse digits from a line using re
# https://stackoverflow.com/questions/4289331/how-to-extract-numbers-from-a-string-in-python
def parseLineDigits(testLine):
    """
    Use regular expressions to extract digits from a string.
    https://stackoverflow.com/questions/4289331/how-to-extract-numbers-from-a-string-in-python

    """
    return re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?", testLine)

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
        One dumpIdy segment, from dumpSegs[], as returned by dumpIdyFileParse()

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

    print(len(dumpSeg))
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
def dumpIdySegsParseX(dumpSegs):
    """
    Extract data from ePS dumpIdy segments into usable form.

    Parameters
    ----------
    dumpSegs : list
        Set of dumpIdy segments, i.e. dumpSegs, as returned by dumpIdyFileParse()

    Returns
    -------
    xr.array
        Xarray data array, containing matrix elements etc.

    int
        Number of blank segments found.


    Example
    -------

    >>> data = parseDumpSegX(dumpSegs)

    """

    dataList = []
    blankSegs = 0

    # Loop over DumpIdy segments, extract data & reformat
    # If blank, skip parser and append blankSegs.
    for dumpSeg in dumpSegs:
        if len(dumpSeg)>4:
            segBlock, attribs = dumpIdySegParse(dumpSeg)
            dataList.append([segBlock[:,0:5].T, segBlock[:,5]+1j*segBlock[:,6], attribs])
            # Switch l,m - with advanced indexing, other methods faster...? https://stackoverflow.com/questions/4857927/swapping-columns-in-a-numpy-array
            # Now done later via pd.MultiIndex.swaplevel()
            # dataList[-1][0][[0,1],:] = dataList[-1][0][[1,0],:]

            # dataList.append([segBlock[:,0:5].T, segBlock[:,5]+1j*segBlock[:,6], attribs])
            # dataList.append([segBlock[:,0:5], segBlock[:,5]+1j*segBlock[:,6], attribs])
        else:
            blankSegs += 1

    # Convert to xarray - ugly loop version, probably a better way to do this!
    #TODO Should:
    #   - be able to loop more cleanly over attribs - set as dict?
    #   - integrate with above loop
    #   - compile into dataSet or dataArray more directly...?
    dataArrays = []
    for data in dataList:
        attribs = data[2]
        QNs = pd.MultiIndex.from_arrays(data[0].astype('int8'), names = attribs[-1][1][0:-1])
        QNs = QNs.swaplevel(0, 1)  # Switch l,m indexes
        #Esyms = pd.MultiIndex.from_arrays([np.array(attribs[0][1]), [attribs[3][1], attribs[5][1]]], names=['E', 'Sym'])
        # pd.MultiIndex.from_tuples([(np.array(attribs[0][1]), [attribs[3][1], attribs[5][1]])], names=['E', 'Sym'])
        Esyms = pd.MultiIndex.from_tuples([(attribs[0][1],attribs[4][1],attribs[5][1],attribs[6][1])],names=[attribs[0][0],attribs[4][0],attribs[5][0],attribs[6][0]])

        #dataArrays.append(xr.DataArray(data[1], coords={'ES': Esyms, 'QN':QNs}, dims = ['ES','QN']))
        # AH - issue is number of labels - can't lable singleton dim it seems, but can expand
        tmp = xr.DataArray(np.asarray(data[1]), coords={'QN':QNs}, dims = ['QN'])
        tmp = tmp.expand_dims({'ES':Esyms})

        # Assign any other attributes - note that some attributes may be dropped when combining arrays below
        for a in attribs:
            tmp.attrs[a[0]] = a[1] # Currently set without units, multiple values here give combine issues below.

        dataArrays.append(tmp)

    # Combine to single xarray
    # Note xarray > v0.12.1
    daOut = xr.combine_nested(dataArrays, concat_dim=['ES'])

    return daOut, blankSegs

# Function for grabbing files or scanning dir for files.
def getFiles(fileIn = None, fileBase = None, fType = '.out'):
    """
    Read ePS file(s) and return results as Xarray data structures.
    File endings specified by fType, default *.out.

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
        File ending for ePS output files, default '.out'


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
        if type(fileIn) is str:
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
        print('\n*** Scanning file(s)')
        print(fList)

    else:
        # Filenames only
        # fList = [f for f in os.listdir(fileBase) if f.endswith(fType)]
        # With full path
        fList = [os.path.join(fileBase, f) for f in os.listdir(fileBase) if f.endswith(fType)]

        # Display message
        print('\n*** Scanning dir')
        print(fileBase)
        print('Found {0} {1} file(s)\n'.format(len(fList), fType))

    return fList

#****** Master function to read a file, or dir, of ePS outputs.

# Some of the logic and methods could do with a revisit/tidy-up here...
#TODO: Add error checking on paths using os.path.isdir/isfile etc.
#TODO: Check/fix paths if incorrectly passed, e.g. https://stackoverflow.com/a/21605790

def readMatEle(fileIn = None, fileBase = None, fType = '.out'):
    """
    Read ePS file(s) and return results as Xarray data structures.
    File endings specified by fType, default *.out.

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
        File ending for ePS output files, default '.out'


    Returns
    -------
    list
        List of Xarray data arrays, containing matrix elements etc. from each file scanned.


    Examples
    --------

    >>> dataSet = readMatEle()  # Scan current dir

    >>> fileIn = r'C:\share\code\ePSproc\python_dev\no2_demo_ePS.out'
    >>> dataSet = readMatEle(fileIn)  # Scan single file

    >>> dataSet = readMatEle(fileBase = r'C:\share\code\ePSproc\python_dev') # Scan dir

    Notes
    -----
    Files are scanned for matrix element output flagged by "DumpIdy" headers.
    Each segment found is parsed for attributes and data (set of matrix elements).
    Matrix elements and attributes are combined and output as an Xarray array.

    """

    # Set files to read, either:
    #   - dir (one only) to scan using os
    #   - file(s) and add to list with full path
    #   - default case is scan current working dir.

    print('*** ePSproc readMatEle(): scanning files for DumpIdy segments (matrix elements)')

    # Call function to check files or scan dir.
    fList = getFiles(fileIn = fileIn, fileBase = fileBase, fType = fType)

    # Loop over fList and scan ePS files
    dataSet = []
    for file in fList:
        print('\n*** Reading ePS output file: ', file)

        # Scan the file and parse segments
        #lines, dumpSegs = dumpIdyFileParse(os.path.join(fileBase, file))
        lines, dumpSegs = dumpIdyFileParse(file)
        data, blankSegs = dumpIdySegsParseX(dumpSegs)

        # Add some additional properties to the output
        fName = os.path.split(file)
        data.name = fName[1]
        data.attrs['file'] = fName[1]
        data.attrs['fileBase'] = fName[0]

        print('Found {0} sets of matrix elements ({1} blank)'.format(len(dumpSegs),blankSegs))

        # Put in a list for now, might want to use Xarray dataset here, and/or combine results from multiple files.
        dataSet.append(data)


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
def readOrb3D(fileIn = None, fileBase = None, fType = '_Orb.dat'):
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
    fList = getFiles(fileIn = fileIn, fileBase = fileBase, fType = fType)

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

    Notes
    -----
    Uses Paulo Herrera's eVTK, see:
    https://pyscience.wordpress.com/2014/09/06/numpy-to-vtk-converting-your-numpy-arrays-to-vtk-arrays-and-files/
    https://bitbucket.org/pauloh/pyevtk/src/default/

    """
    fOut = []
    for file in dataSet:
        # Set grid, convert to Cart if necessary, assuming that grid won't be larger than 10 Angs
        if (len(file[2][0]) > 50):
            # Convert to Cart grid for plotting
            # TODO: Investigate use of sph grid here - should be cleaner.
            # TODO: Investigate recreating mesh in Paraview, rather than saving to file.
            [T,R,P] = np.meshgrid(file[2][1], file[2][0], file[2][2])
            T = (T*np.pi/180) #-np.pi/2
            P = P*np.pi/180
            x = R*np.sin(P)*np.cos(T)
            z = R*np.cos(P)
            y = R*np.sin(P)*np.sin(T)
        else:
            [x,y,z] = np.meshgrid(file[2][1], file[2][0], file[2][2])

        # Save single dataset
        # Function info: https://bitbucket.org/pauloh/pyevtk/src/default/src/hl.py
        # gridToVTK("./ePStest3", x, y, z, pointData = {'data': data[0]})

        # Save full dataset
        # TODO: Check number of segments and save multiple
        # segs = round(len(file[3])/3)
        fOut.append(gridToVTK(file[0][:-4], x, y, z, pointData = {'Re': file[3][0], 'Im': file[3][1], 'Abs': file[3][2]}))

    print("{0} files written to vtk format.".format(len(fOut)))

    return fOut