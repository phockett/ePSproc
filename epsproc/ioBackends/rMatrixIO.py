"""
ePSproc IO routines for R-Matrix codes

- Read R-Matrix matrix elements
    - Use existing IO style (seg passing style), although may be worth rewriting for speed at some point.
- Convert format to match ePolyScat conventions (spherical basis, complex spherical harmonics).
    - Main sph routines (sphCalc) now updated with real harmonics.
    - Util funcs added to sphFuncs.sphConv.

v1  20/09/22

"""


import xarray as xr
import pandas as pd
import numpy as np

from pathlib import Path
from io import StringIO
import re

# from epsproc.IO import fileParse   # 03/10/22 - circular import errors currently, should reorg base IO to fix.
                                     # For now use lambda func work-around, see https://stackoverflow.com/a/42114399
import epsproc.IO
# fileParse = lambda: epsproc.IO.fileParse
IO = lambda: epsproc.IO
# fileParse = IO().fileParse

from epsproc.util.conversion import conv_ev_atm



#******* LOW-LEVEL IO

def RmatChannelSeg(fileName, verbose = 1):
# Simple wrapper for general fileParse function, ePS dumpIdy segments
# ADAPTED FROM def dumpIdyFileParse(fileName, verbose = 1):
    """
    Parse an R-matrix file for channel list segments.

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

    startPhrase = "Channel    Targ."  # startPhrase skips whitespace
    endPhrase = "  Channel     Electron Energy"   # Note endPhrase keeps leading whitespace.

    (lines, dumpSegs) = IO().fileParse(fileName, startPhrase, endPhrase, verbose = verbose) # , '>')

    # NOTE - with current code dumpSegs has one final blank segment
    if verbose:
        print('Found {0} Channel list segments.'.format(len(dumpSegs) - 1))

    return lines, dumpSegs[:-1]


def RmatHeader(fileName, lines):
    """
    Read header from R-matrix dipoles file.

    Parameters
    ----------
    fileName : str
        File to read (file in working dir, or full path)

    lines : int, tuple, list or set
        Lines to read. If int read 0:lines, if tuple read range(line[0],line[1]), if list read line[0]:line[1]
        Note these are assumed to be 1-indexed to match file line numbers.

    verbose : bool or int, optional
        If true, print segment info.


    Returns
    -------

    list
        linesOut, [[i, line]] list of lines read from file.


    """

    # Input handling.
    # Note use of set for speed, inspired by https://stackoverflow.com/a/2082131
    if type(lines) is int:
        lines = set(range(1,lines+1))  # +1 to include final line # in list.

    elif type(lines) is tuple:
        lines = set(range(lines[0],lines[1]+1))      # +1 to include final line # in list.

    elif type(lines) is list:
        lines = set(lines)

    with open(fileName,'r') as f:
        linesOut = [[i+1,x] for i, x in enumerate(f) if i+1 in lines]  # i+1 to use 1-indexed line numbers (0-indexed from enumerate())

    return linesOut


def RmatDipoles(fileName, verbose = 1):
# Simple wrapper for general fileParse function, ePS dumpIdy segments
# ADAPTED FROM def dumpIdyFileParse(fileName, verbose = 1):
    """
    Parse an R-matrix file for channel list segments.

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

#     startPhrase = "Channel     Electron Energy"
    startPhrase = "Channel     Electron Energy       real Y"
    endPhrase = None

    (lines, dumpSegs) = IO().fileParse(fileName, startPhrase, endPhrase, verbose = verbose) # , '>')

    # NOTE - with current code dumpSegs has one final blank segment EXCEPT FOR NONE type ending!
    if verbose:
        print('Found {0} dipole segments (sets of matrix elements).'.format(len(dumpSegs)))

    return lines, dumpSegs


#****** Parsers

# Parse a DumpIdy segment
#TODO: More attribs, and convert attribs to dict, or other more structured format.
#TODO: Error checking for missing matrix elements or IO issues - see Matlab code.
#TODO: Convert attribs to dict for simpler assignments later
def RmatDipolesSegParse(dumpSeg):
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
    for testLine in dumpSeg[2:-1]:    #dumpSeg[4000:4020]:    # dumpSeg[2:-1]:

#        vals = parseLineDigits(testLine[2])
#        lineOut = []
#        [lineOut.append(float(val)) for val in vals] # Convert to float
#        rawIdy.append(lineOut)

        # Explicit test for empty lines here.
        # In test R-matrix output files these appear between different channels.
        # Prevents NP array stacking if not skipped.
        if not testLine[2].isspace():
            tmp=np.genfromtxt(StringIO(testLine[2].replace('D','E')))  # Quick replace to handle D notation - might be a cleaner way?
                                                                       # Possibly reading whole file - headers, then lambda fn, e.g. https://stackoverflow.com/a/65945191
                                                                       # Likely will be quite a bit quicker.
    #         if tmp.any():
            rawIdy.append(tmp)

    # Parse header lines (structured, selected by line #)
#     attribs.append(['E', np.float(parseLineDigits(dumpSeg[3][2])[0]), 'eV'])
#     attribs.append(['Ehv', np.float(parseLineDigits(dumpSeg[4][2])[0]), 'eV'])
#     SF = np.genfromtxt(parseLineDigits(dumpSeg[5][2]))
#     SF = SF[0] + SF[1]*1j
#     attribs.append(['SF', SF, 'sqrt(MB)'])
#     attribs.append(['Lmax', np.int(parseLineDigits(dumpSeg[10][2])[0]), ''])

    # Parse symmetries - multiple .spilt(), or use re?
    # attribs.append(['Syms', parseLineDigits(dumpSeg[7][2]), ''])
#     symList = dumpSeg[7][2].split()
    # Append as nested lists
#    sym = []
#    [sym.append([symList[n], symList[n+1].strip('=')]) for n in range(1,6,2)]
#    attribs.append(['Syms', sym, ''])
    # Append as indendent lists
#     [attribs.append([symList[n], symList[n+1].strip('=')]) for n in range(1,6,2)]

#     attribs.append(['QNs', dumpSeg[11][2].split(), ''])

    return np.asarray(rawIdy), attribs



def RmatDipolesSegParseX(dumpSegs, channelSegs, headerSegs, fileName = None):
    dipoles = dumpSegs
    # Complex form - not sure if col inds are guarenteed here?
    dipolesC = dipoles[:, 2::2] + 1j*dipoles[:, 3::2]

    # Test with Channel as main index
    # xrTest = xr.DataArray(dipolesC, coords={'Channel':dipoles[:,0].astype('int8'), 'Pol':['Y','Z','X'], 'E':('Channel',dipoles[:,1])}, dims = ['Channel','Pol'])

    # Test with Multiindex - OK
    # Seems to be simplest thing to do, since it sorts out the restacking.

    ind = pd.MultiIndex.from_arrays([dipoles[:,0].astype('int8'),dipoles[:,1]], names = ['Channel','E'])
    xrDipoles = xr.DataArray(dipolesC, coords={'CE':ind, 'Pol':['Y','Z','X']}, dims = ['CE','Pol']).unstack()

    # xrTest  #.unstack()


    # With E, then reshape & set channel, ugh.
    # xrTest = xr.DataArray(dipolesC, coords={'Pol':['Y','Z','X'], 'E':dipoles[:,1]}, dims = ['E','Pol'])

    # xrTest.reshape(xrTest.E.size)
    # 'Channel':dipoles[:,0].astype('int8'),

    # xrTest.unstack()
    # xrTest.expand_dims('E')


    #**** Add channel info...

    # Split digits to list with parser - OK
    # [ep.IO.parseLineDigits(testLine[2].replace('D','E')) for testLine in channelSegs[0][2:]]

    # Split digits to list with whitespace - OK
    # [testLine[2].replace('D','E').split() for testLine in channelSegs[0][2:]]

    # Convert to NP by line - OK
    # np.asarray([np.genfromtxt(StringIO(testLine[2].replace('D','E'))) for testLine in channelSegs[0][2:]])

    # Convert to PD
    # pd.DataFrame([ep.IO.parseLineDigits(testLine[2].replace('D','E')) for testLine in channelSegs[0][2:]])  # OK

    # pd.DataFrame([testLine[2].replace('D','E').split() for testLine in channelSegs[0][2:]])  # OK

    # With col names - ASSUMES FIXED FORMAT!
    # NOTE need to set dtypes here - objects at the moment.
    cHeader = channelSegs[0][0][2].split()
    cNames = cHeader[0:-2] + [cHeader[-1]]

    channelDF = pd.DataFrame([testLine[2].replace('D','E').split() for testLine in channelSegs[0][2:]], columns = cNames).dropna()
    # channelDF = pd.DataFrame([testLine[2].replace('D','E').split() for testLine in channelSegs[0][2:]], columns = cNames).dropna()
    # .set_index('Channel')  # OK

    # Tidy up formatting
    # channelDF = channelDF.iloc[:,0:-1].astype(int)
    # channelDF.to_numeric(df['DataFrame Column'],errors='coerce')

    # Ugh, explictly set dtypes for conversion - alternative is to loop over cols with pd.to_numeric(), but this just casts all to float
    # channelDFcast = channelDF.apply(pd.to_numeric, axis=1)  # All to floats
    cDtypes = {k:int for k in cHeader[0:-2]}  # Set dict of cols:dtype
    cDtypes.update({cHeader[-1]:float})

    # Assign dtypes & reindex
    # channelDF = channelDF.astype(dtype=cDtypes).set_index(['Channel', 'l', 'm'])   # Include Channel?
    channelDF = channelDF.astype(dtype=cDtypes).set_index(['l', 'm'])  # Skip channel

    # Test XR with new multindex...
    # Config as before
    ind = pd.MultiIndex.from_arrays([dipoles[:,0].astype('int8'),dipoles[:,1]], names = ['Cind','Eh'])
    xrDipoles = xr.DataArray(dipolesC, coords={'CE':ind, 'Pol':['Y','Z','X']}, dims = ['CE','Pol']).unstack()

    # Update Cind to new multiindex
    # May also use something like (for items instead of existing index): cInd = pd.MultiIndex.from_frame(channelDF)

    # Case with Channel included - not sure if required?
    # xrDipoles.coords['Cind'] = channelDF.index
    # xrDipoles = xrDipoles.unstack()  #.stack({'LM':['l','m']})
    # Should now be able to drop Channel index?  xrDipoles.drop('Channel') ?
    # But having some issues with this as index coord.
    # Bit more tidy up testing - should check existing code techniques here!
    # xrDipoles.drop('Channel').squeeze(drop=True)
    # xrDipoles.drop_dims('Channel')
    # xrDipoles.reset_coords('Channel')   # This fails as Channel is index, also note LM=28, CHannel=16, so needs cleaning up!
    # xrTidy = xrDipoles.dropna(how='all', dim='LM').drop('Channel') # Might be OK, although still shows Channel.  # .reset_coords('Channel')    #.reset_coords('Channel',drop=True)


    # For l,m index only
    xrDipoles.coords['Cind'] = channelDF.index
    xrDipoles = xrDipoles.rename({'Cind':'LM'})


    # E conversion H > eV
    # cf. Pint with Xarray - may want to implement, https://pint-xarray.readthedocs.io/en/stable/creation.html
    # Also basic unit labels with XR: https://docs.xarray.dev/en/stable/getting-started-guide/quick-overview.html#attributes
    xrDipoles = xrDipoles.assign_coords(Eke=('Eh', conv_ev_atm(xrDipoles.Eh).values)).swap_dims({'Eh':'Eke'})
    xrDipoles['Eh'].attrs['units'] = 'au'
    xrDipoles['Eke'].attrs['units'] = 'eV'

    # Additional attrs - should add some further processing here
    # xrDipoles.attrs['header'] = headerSegs
    xrDipoles.attrs['header'] = ' '.join([str(x[1]) for x in headerSegs])
    xrDipoles.attrs['file'] = fileName

    # Add pol geom(?) from file name.
    # Not sure if this is general, so use try/except for lazy handling.
    if fileName is not None:
        suff = Path(fileName).suffix
        matches = re.search(f'(\w){suff}$',fileName)

        if matches:
            try:
                xrDipoles = xrDipoles.expand_dims({'LF':[matches.group(1)]})
            except Exception as e:
                print(f"*** Warning: Failed to expand dims with filename polarization '{matches.group(1)}', caught {type(e).__name__}: {e}")

    # TODO
    # - Additional phase terms - see 'Coulomb phase' section below
    # - Convert E values & set with .assign_coords(Eke=('E',Eke))
    # - Thresholding?
    # - Real > complex harmonics.
    # - Additional indexes? Pol cart > spherical?

    return xrDipoles

#******************* SPH ROUTINES NOW IN sphFuncs.sphConv module.
#
# import numpy as np
# import xarray as xr
#
# from epsproc.util.listFuncs import genLM
#
# try:
#     import pyshtools as pysh
# except ImportError as e:
#     if e.msg != "No module named 'pyshtools'":
#         raise
#     print('* pyshtools not found, SHtools functions not available. If required, run "pip install pyshtools" or "conda install -c conda-forge pyshtools" to install.')
#
#
# # Set SHtools object and init with set_coeffs
# # x.set_coeffs(values, ls, ms)
# # NOTE THIS CASTS TO REAL for real harmonic, so not sure how to convert dipoles here - will only apply to real part...?
#
# # npTab = dfLong.to_numpy()
# # clm = pysh.SHCoeffs.from_zeros(lmax = 6)   # Defaults to kind = 'real'
#
# def SHcoeffsFromXR(dataIn, kind = 'real'):
#     """
#     Xarray Spherical Harmonic coeffs to SHtools
#
#     MAY HAVE BETTER VERSION ELSEWHERE?
#
#     This will only work for 1D array -- will need to slice/group for more general handling.
#
#     TODO: update with general l,m dim names and handling, see checkSphDims()
#     """
#
#     if ('kind' in dataIn.attrs.keys()):
#         kind = dataIn.attrs['kind']
#
#     if ('harmonics' in dataIn.attrs.keys()):
#         kind = dataIn.attrs['harmonics']['kind']
#
# #     else:
# #         kind = kind
#
#     clm = pysh.SHCoeffs.from_zeros(lmax = dataIn.l.max().values, kind = kind)
#     clm.set_coeffs(dataIn.values, dataIn.l.astype(int), dataIn.m.astype(int))  # NEEDS (values, ls, ms)
#
#     return clm
#
#
#
#
# def XRcoeffsFromSH(dataIn):
#     """
#     Convert SHtools coeffs object to Xarray
#
#     NOTE: see ep.genLM and relate utils.
#     ep.sphCalc.setBLMs  is similar too, although from list not SHtools object.
#
#     """
#
#     # Set from inds
# #     # Generate indexes
# #     lmax = dataIn.lmax
# #     LMlist = genLM(lmax).T
# #     ind = pd.MultiIndex.from_arrays(LMlist, names = ['l','m'])
#
# #     # Set +ve and -ve arrays to match SHtools convention [+/-,l,m] index
# #     posMind = LMlist[1,:] >-1
# #     negMind = LMlist[1,:] <1
# #     indPos = pd.MultiIndex.from_arrays(LMlist[:,posMind], names = ['l','m'])
# #     indNeg = pd.MultiIndex.from_arrays(LMlist[:,negMind], names = ['l','m'])
#
#
# #     return indPos
#
#
# #     # Direct to XR with sign
# #     import numpy as np
# #     lmXRTest = xr.DataArray(clmC.coeffs, coords={'sign':[0,1],'l':np.arange(0,clmC.coeffs.shape[1]), 'm':np.arange(0,clmC.coeffs.shape[2])})   #, coords={'l':})
# #     lmXRTest
#
#
#     # Split and restack - might be easier?
#     lmXRTest1 = xr.DataArray(dataIn.coeffs[0,:,:], coords={'l':np.arange(0,dataIn.coeffs.shape[1]), 'm':np.arange(0,dataIn.coeffs.shape[2])})   #, coords={'l':})
#     lmXRTest2 = xr.DataArray(dataIn.coeffs[1,:,:], coords={'l':np.arange(0,dataIn.coeffs.shape[1]), 'm': -1*np.arange(0,dataIn.coeffs.shape[2])})   #, coords={'l':})
#
#     # Issues with duplicate l,m=0,0 case here
#     # lmXRTest = xr.concat([lmXRTest1,lmXRTest2], dim='m', compat='no_conflicts', join='override')
#
#     # Gives NaNs for duplicate case
#     lmXRTest = xr.concat([lmXRTest1,lmXRTest2.where(lmXRTest2.l>0).where(lmXRTest2.m<0)], dim='m')  #, compat='no_conflicts', join='override')
#
#     # Clean up - HAVE CODE FOR THIS ELSEWHERE?
#     lmXRTestClean = lmXRTest.where(np.abs(lmXRTest.m)<=lmXRTest.l, drop=True)
#     lmXRTestClean = lmXRTestClean.stack({'LM':['l','m']}).dropna(dim='LM')
#
#     # Threshold?
#     lmXRTestClean = lmXRTestClean.where(np.abs(lmXRTestClean)>1e-4).dropna(dim='LM')
#
#     lmXRTestClean.attrs['SH'] = dataIn
#     lmXRTestClean.attrs['kind'] = dataIn.kind
#
#     return lmXRTestClean
#
#
# def checkSphDims(dataIn, keyDims = {'LM':['l','m']}):
#     """
#     Very basic dim handling/checking for spherical harmonic Xarrays.
#
#     Parameters
#     ----------
#     dataIn : Xarray
#         Xarray containing values to convert.
#         Must contain keyDims (stacked or unstacked).
#         TODO: generalise this, should get from :py:func:`epsproc.util.listFuncs.getRefDims()`.
#
#     keyDims : dict, optional, default = {'LM':['l','m']}
#         Defines dim names for harmonic rank and order, and stacked dim name.
#
#
#     Returns
#     -------
#     dataIn : Xarray
#         As dataIn, but with additional .attrs['harmonics'] defining dim relations, and unstacked along keyDims[stackDim].
#
#
#     TODO: should use ep.util.misc.checkDims here?
#     TODO: check other dim handling functionality, may be reinventing things here.
#
#     """
#
#     LMStackFlag = False
#     stackDim = list(keyDims.keys())[0]
#     dimList = keyDims[stackDim]
#     lDim = keyDims[stackDim][0]
#     mDim = keyDims[stackDim][1]
#
#     if stackDim in dataIn.dims:
#         dataIn =  dataIn.unstack(stackDim)
#         LMStackFlag = True
#
#     if 'harmonics' not in dataIn.attrs.keys():
#         dataIn.attrs['harmonics'] = {}
#
#     dataIn.attrs['harmonics'].update({k:v for k,v in locals().items() if k != 'dataIn'})
#
#     return dataIn
#
#
#
# def sphRealConvert(dataIn, method = 'std', keyDims = {'LM':['l','m']}, incConj = True, rotPhase = None):
#     """
#     Convert real harmonics to complex form.
#
#     Xarray for input and outputs.
#
#     See also `PEMtk.sym.symHarm.setCoeffsSH`
#
#     Parameters
#     ----------
#     dataIn : Xarray
#         Xarray containing values to convert.
#         Must contain keyDims (stacked or unstacked).
#         TODO: generalise this, should get from :py:func:`epsproc.util.listFuncs.getRefDims()`.
#
#     method : str, optional, default = 'std'
#         - 'std': Set +/-m terms from |m| real harmonics as per wikipedia defn, https://en.wikipedia.org/wiki/Spherical_harmonics#Real_form
#                 Note: this form currently assumes only +M or -M for any given term, and also sets conjugates upon conversion if incConj=True.
#                 This may be phase-rotated relative to SHtools convention.
#         - 'sh': Set terms as per SHtools definitions, https://shtools.github.io/SHTOOLS/complex-spherical-harmonics.html
#                 Note this form has purely real or purely imaginary terms.
#
#     keyDims : dict, optional, default = {'LM':['l','m']}
#         Defines dim names for harmonic rank and order, and stacked dim name.
#
#     incConj : bool, optional, default = True
#         Set opposite sign M term as conjugate in conversion.
#         Generally should be True, but may cause issues in cases where both +/-M are defined in input.
#
#     rotPhase : float, optional, default = None
#         If not None, apply additional phase rotation by exp(-i*rotPhase), i.e. rotate phase about Z-axis.
#         To match SHtools defn. set rotPhase = -np.pi/2 which will flip (x,y) axes.
#
#     Returns
#     -------
#
#     dataC : Xarray
#         Results of conversion.
#         All parameters set in dataC['harmonics']
#
#
#     TODO
#
#     - Generalise rotations, currently only supports theta (about Z-axis) to allow matching to SHtools conventions.
#         See ep.sphCalc.TKQarrayRotX, which does work but needs some tidying up.
#         E.g.    # Test existing frame rot... quite old code
#                 testBLMs = sphRealConvert(lmReal, method='std').rename({'LM':'BLM'})
#                 testBLMs.attrs['dataType'] = 'BLM'
#                 BLMrot, *_ = ep.TKQarrayRotX(testBLMs, ep.setPolGeoms())
#
#     """
#
#     dataCalc = dataIn.copy()
#
#     #*** Basic dim handling
#     # TODO: should use checkDims here
#     # LMStackFlag = False
#     # stackDim = list(keyDims.keys())[0]
#     # lDim = keyDims[stackDim][0]
#     # mDim = keyDims[stackDim][1]
#     #
#     # if stackDim in dataCalc.dims:
#     #     dataCalc =  dataCalc.unstack(stackDim)
#     #     LMStackFlag = True
#
#     # Functional version - returns above params to dataCalc.attrs['harmonics']
#     dataCalc = checkSphDims(dataCalc, keyDims)
#     mDim = dataCalc.attrs['harmonics']['mDim']  # For brevity below.
#
#     #*** METHODS
#     if (method == 'std-v1') or (method == 'std'):
#         # Set +/-m terms from |m| real harmonics as per
#         # https://en.wikipedia.org/wiki/Spherical_harmonics#Real_form
#
#         # Compute for all m, subselect +/- later.
#         # plusMcase = ((-1)**np.abs(dataCalc.m))/np.sqrt(2)*(dataCalc + 1j*dataCalc)
#         # plusMcase = ((-1)**np.abs(dataIn.m))/np.sqrt(2)*(dataIn + 1j*dataIn)  # For -ve m may need abs(m) here
#         negMcase = (1/np.sqrt(2))*(dataCalc - 1j*dataCalc)
#         plusMcase = ((-1)**np.abs(dataCalc[mDim]))/np.sqrt(2)*(dataCalc + 1j*dataCalc)
#
#         # Concat from various parts - OK.
#         nonZeroM = xr.concat([plusMcase.where(plusMcase[mDim]>0, drop=True),  negMcase.where(negMcase[mDim]<0, drop=True)], dim=mDim)  # Keep all m in both cases
#
#         # Set conjugate pairs?
#         if incConj:
#             dataC = xr.concat([dataCalc.sel({mDim:0}), nonZeroM, sphMswap(nonZeroM)], dim=mDim)
#         else:
#             dataC = xr.concat([dataCalc.sel({mDim:0}), nonZeroM], dim=mDim)
#
#     elif method == 'std-v2':
#         # Set by coord assignment... have some norm issues with v1?
#         # Test with subselection & different stacking...
#         plusMterms = dataCalc.where(dataCalc[mDim]>0)
#         negMterms = dataCalc.where(dataCalc[mDim]<0)
#
#         plusMcase = ((-1)**np.abs(plusMterms[mDim]))/np.sqrt(2)*(plusMterms + 1j*plusMterms)
#         # plusMcase = ((-1)**np.abs(dataIn.m))/np.sqrt(2)*(dataIn + 1j*dataIn)  # For -ve m may need abs(m) here
#         negMcase = (1/np.sqrt(2))*(negMterms - 1j*negMterms)
#
#         # Concat from various parts - OK.
#         nonZeroM = xr.concat([plusMcase.where(plusMcase[mDim]>0, drop=True),  negMcase.where(negMcase[mDim]<0, drop=True)], dim=mDim)  # Keep all m in both cases
#
#         # Set conjugate pairs?
#         if incConj:
#             dataC = xr.concat([dataCalc.sel({mDim:0}), nonZeroM, sphMswap(nonZeroM)], dim=mDim)
#         else:
#             dataC = xr.concat([dataCalc.sel({mDim:0}), nonZeroM], dim=mDim)
#
#
#     elif method == 'sh':
#         # Set by coord assignment... Should match SHtools case.
#         dataC = xr.zeros_like(dataCalc)
#         # dataC = dataC.where(dataC[mDim]>-1, 1/np.sqrt(2)*(- 1j*dataCalc))   # -m case
#         dataC = dataC.where(dataC[mDim]>-1,((-1)**np.abs(dataCalc[mDim]))/np.sqrt(2)*(1j*dataCalc))   # -m case
#         # dataC3 = dataC3.where(dataC3.m<1,((-1)**np.abs(dataIn.m))/np.sqrt(2)*(dataIn + 1j*dataIn))   # +m case
#         dataC = dataC.where(dataC[mDim]<1,1/np.sqrt(2)*(dataCalc))   # +m case
#         dataC = dataC.where(dataC[mDim]!=0,dataCalc)  # m=0 case
#
#         if incConj:
#             # Generate +/-m pairs from above
#             dataC = xr.concat([dataC, sphMswap(dataC.where(dataC[mDim]!=0))], dim=mDim)  #, coords='minimal')
#             # dataCp = dataCp.sortby(keyDims[stackDim]).stack(keyDims).dropna(dim=stackDim,how='all')
#
#     else:
#         print(f"*** sphRealConvert(method = {method}) not recognised.")
#
#         return 0
#
#
#     #*** Tidy up
#     if rotPhase is not None:
#         # To match SHtools defn. set rotPhase = -np.pi/2 which will flip (x,y) axes.
#         dataC = dataC * np.exp(-1j*rotPhase)
#
#     # dataC = dataC.sortby(keyDims[stackDim]).stack(keyDims).dropna(dim=stackDim,how='all')  # Sort dims, stack and clean up.
#     dataC = dataC.sortby(dataCalc.attrs['harmonics']['dimList']).stack(keyDims).dropna(dim=dataCalc.attrs['harmonics']['stackDim'],how='all')  # Sort dims, stack and clean up.
#
#     dataC.attrs['harmonics'] = {'dtype':'Complex harmonics',
#                                 'kind':'complex',
#                                 'method': method,
#                                 'incConj':incConj,
#                                 'keyDims':keyDims}
#
#     return dataC
#
#
# # DEV CODE
# # def sphRealConvert(dataIn):
# #     """
# #     Convert real harmonics to complex form.
# #
# #     Xarray for input and outputs.
# #
# #     See also `PEMtk.sym.symHarm.setCoeffsSH`
# #
# #     """
# #
# #     # Set +/-m terms from |m| real harmonics as per
# #     # https://en.wikipedia.org/wiki/Spherical_harmonics#Real_form
# #
# #     # Basic test - assume |m| > +/-m terms.
# #     # May be incorrect phase for -m real harmonics?
# #     # Note need abs(m), otherwise (-1)**m may throw errors
# #     plusMcase = ((-1)**np.abs(dataIn.m))/np.sqrt(2)*(dataIn + 1j*dataIn)
# #     # plusMcase = ((-1)**np.abs(dataIn.m))/np.sqrt(2)*(dataIn + 1j*dataIn)  # For -ve m may need abs(m) here
# #     negMcase = (1/np.sqrt(2))*(dataIn - 1j*dataIn)
# #     negMcase = negMcase.unstack('LM')
# #     negMcase.coords['m'] = -negMcase.m
# #     # negMcase = negMcase.stack({'LM':['l','m']})
# #
# #     # For case with ONLY + or -m in real harmonics (should be general...?) need to set by symmetry relations
# #     # TODO: is this general? Should be since it's just a phase convention?
# #     # THIS MATCHES WIKI defn., but in testing gives conjugate answer - might be accidental switch in numerics?
# #     minM = dataIn.m.min()
# #     if minM <0:
# #         negMcase2 = (1/np.sqrt(2))*(dataIn - 1j*dataIn)
# #         plusMcase2 = np.conj(((-1)**np.abs(dataIn.m)) *negMcase2)
# #
# #         negMcase2 = negMcase2.unstack('LM')
# #         plusMcase2 = plusMcase2.unstack('LM')
# #
# #         plusMcase2.coords['m'] = -plusMcase2.m
# #
# #
# #     if minM==0:
# #         plusMcase2 = ((-1)**dataIn.m)/np.sqrt(2)*(dataIn + 1j*dataIn)
# #
# #         negMcase2 = np.conj(((-1)**np.abs(dataIn.m)) *plusMcase2)
# #
# #         plusMcase2 = plusMcase2.unstack('LM')
# #         negMcase2 = negMcase2.unstack('LM')
# #         negMcase2.coords['m'] = -negMcase2.m
# #
# #
# #     # Concat - BASICALLY WORKING, just need to drop/fix m=0 terms.
# #     # UPDATE - now fixed with .where()
# #
# #     # dataC = plusMcase + negMcase   # Keeps m=0 terms only!
# #
# #     # Concat from various parts - OK.
# #     dataC = xr.concat([dataIn.unstack('LM').sel(m=0), plusMcase.where(plusMcase.m>0, drop=True).unstack('LM'),
# #                        negMcase.where(negMcase.m<0, drop=True)], dim='m')
# #     #Tidy up
# #     dataC = dataC.stack({'LM':['l','m']}).dropna(dim='LM',how='all')
# #     dataC.attrs['dtype'] = 'Complex harmonics'
# #     dataC.attrs['kind'] = 'complex'
# #     # dataC = dataC.stack({'LM':['l','m']}).dropna(dim='LM',how='all')
# #
# #
# # #     return dataC, dataIn, plusMcase2, negMcase2
# #
# #     # Test 2nd form...
# #     # Concat from various parts
# #
# #     # V1 - kills some terms due to removal of +/-m coords? May be able to fix with coords = something?
# # #     dataC2 = xr.concat([dataIn.unstack('LM').sel(m=0), plusMcase2.where(plusMcase2.m>0, drop=True),
# # #                        negMcase2.where(negMcase2.m<0, drop=True)], dim='m')  #, coords='minimal')
# #
# #     # V2 - correctly treats m=0 and keeps all other terms, not sure if other terms ARE CORRECT HOWEVER
# #     dataC2 = xr.concat([dataIn.unstack('LM').sel(m=0),
# #                     plusMcase2.where(plusMcase2.m != 0),
# #                     negMcase2.where(negMcase2.m != 0)], dim='m')  #, coords='minimal')
# #
# #     #Tidy up
# #     dataC2 = dataC2.stack({'LM':['l','m']}).dropna(dim='LM',how='all')
# #     dataC2.attrs['dtype'] = 'Complex harmonics'
# #     dataC2.attrs['kind'] = 'complex'
# #
# #     return dataC, dataC2
# #
# #     # V3 - swap on coords only
# #     dataC3 = xr.zeros_like(dataIn)
# #     dataC3.where(dataC3.m<0,(1/np.sqrt(2))*(dataIn - 1j*dataIn))
# #
# #
# #
# # #     # V4 - defn. from https://shtools.github.io/SHTOOLS/complex-spherical-harmonics.html
# # #     # This defines coeff transforms, so should be correct.
# # #     plusMcase = 1/np.sqrt(2)*(dataIn - 1j*dataIn)
# # #     # plusMcase = ((-1)**np.abs(dataIn.m))/np.sqrt(2)*(dataIn + 1j*dataIn)  # For -ve m may need abs(m) here
# # #     negMcase = (1/np.sqrt(2))*(dataIn - 1j*dataIn)
# # #     negMcase = negMcase.unstack('LM')
# # #     negMcase.coords['m'] = -negMcase.m
# #
#
#
#
# def sphMswap(dataIn, cleanOutput = True):
#     """
#     Compute +/-m sign switch as conjugate harmonics from Xarray.
#
#     Ylm = conj((-1)^m * Yl,-m)
#
#     Basic dim handling with checkSphDims(dataIn) implemented, but should generalise.
#
#     """
#
#     dataCalc = dataIn.copy()
#
#     # LMStackFlag = False
#     # if 'LM' in dataCalc.dims:
#     #     dataCalc =  dataCalc.unstack('LM')
#     #     LMStackFlag = True
#
#     # Functional version - returns above params to dataCalc.attrs['harmonics']
#     # Only run if not already set (e.g. already run in sphRealConvert(), including keyDims)
#     # if (not 'harmonics' in dataCalc.attrs.keys()) or (not 'mDim' in dataCalc.attrs['harmonics'].keys()):
#     #     dataCalc = checkSphDims(dataCalc)   # Just run with default keyDims - needs rationalisation.
#
#     # Always run? This avoids issues with passed attrs in some cases (e.g. if set on a different dataIn style)
#     # TODO: overwrite option?
#     dataCalc = checkSphDims(dataCalc)   # Just run with default keyDims - needs rationalisation.
#
#     mDim = dataCalc.attrs['harmonics']['mDim']  # For brevity below.
#
#     dataOut = np.conj((-1)**np.abs(dataCalc[mDim]) * dataCalc)    # Conj resultant
#     # dataOut = (-1)**np.abs(dataCalc.m) * np.conj(dataCalc)  # Conj inputs only
#     dataOut[mDim] = -1 * dataOut[mDim]
#
#     if cleanOutput:
#         dataOut = dataOut.dropna(dim=mDim,how='all')
#
#     if dataCalc.attrs['harmonics']['LMStackFlag']:
#         dataOut = dataOut.stack(dataCalc.attrs['harmonics']['keyDims'])
#
#         # Tidy again?
#         if cleanOutput:
#             dataOut = dataOut.dropna(dim=dataCalc.attrs['harmonics']['stackDim'],how='all')
#
#     return dataOut
