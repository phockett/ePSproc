"""
ePSproc classes IO function wrappers

16/10/20    Methods for/from base class.
"""

from pathlib import Path
import pprint

# Local functions
from epsproc import readMatEle, headerFileParse, molInfoParse, multiDimXrToPD, plotTypeSelector, matEleSelector
from epsproc.util.summary import getOrbInfo, molPlot
from epsproc.util.env import isnotebook

def scanFiles(self, dataPath = None, fileIn = None, reset = False, keyType = 'orb', **kwargs):
    """
    Scan ePS output files from a dir for multiple data types. Sort data, and set to list/dict/Xarray structures.

    Adapted from https://phockett.github.io/ePSdata/XeF2-preliminary/XeF2_multi-orb_comparisons_270320-dist.html

    Current implementation:
    - Read XS and matrix elements from source files, sort to Xarrays (one per file and data type), uses uses :py:func:`epsproc.readMatEle`.
    - Stack by Eke for multi-file E-chunked jobs.
    - Read additional data for jobs (uses :py:func:`epsproc.headerFileParse` and :py:func:`epsproc.molInfoParse`).
    - Sort data to lists by data type, and dict with keys per file/job (self.data).
    - Dict should be final data to use (note - can't get heterogenous data types & dims to work well for Xarray Dataset, but this may change.)

    TODO:
    - convert outputs to Xarray dataset. Did this before, but currently missing file (on AntonJr)! CHECK BACKUPS - NOPE.
    - Confirm HV scaling - may be better to redo this, rather than correct existing values?

    - Fix xr.dataset: currently aligns data, so will set much to Nan if, e.g., different symmetries etc.
    Change to structure as ds('XS','matE') per orb, rather than ds('XS') and ds('matE') for all orbs?
    This should also be in line with hypothetical base dataclass, which will be per orb by defn.

    Parameters
    -----------
    dataPath : str or Path object, optional, default = None
        Set dir to scan.
        Default is to use self.job['fileBase'] as set at init.

    reset : bool, optional, default = False
        If False, new data will be appended to any existing data.
        If True, existing data will be removed.
        This allows for persistence over multiple calls, e.g. reading multiple dirs.

    keyType : str, optional, default = 'orb'
        'orb': Use orbital labels as dataset keys
        'int': Use integer labels as dataset keys (will be ordered by file read)
        Any other setting will result in key = keyType, which can be used to explicitly pass a key (e.g. in multijob wrapper case). This should be tidied up.

    **kwags : optional
        Additional args passed to base ep.readMatEle() routine.
        For non-standard files, set as required by readMatEle(), and some post-processing will be skipped.
        E.g. to read R-matrix, run with `self.scanFiles(fType=fType, recordType='Rmat', stackDim = 'LF', keyType='name')`


    07/10/22    Added **kwargs and alternative IO post-processing for other file types.
                Quite messy, should tidy main loop and integrate.

    """

    # Allow dir & file override here for subclasses/independent use
    if fileIn is None:
        fileIn = self.job['fileIn']  # This will be set to None if not set, which is fine for readMatEle

    # Set to None if a filelist is provided
    if (dataPath is None) and (fileIn is None):
        dataPath = self.job['fileBase']
#         else:
#             self.job['fileBase']


    # Remove existing data?
    if reset:
        self.data = {}
        self.jobNotes = []

    # For Xarray Dataset stacking...
    # dsXS = xr.Dataset()  # Set blank dataset
    # dsMatE = xr.Dataset()

    # 07/10/22 - added general file wrapper here for non-default cases
    if 'recordType' in kwargs.keys():
        dataSetMatE = readMatEle(fileBase = dataPath, fileIn = fileIn, verbose = self.verbose['sub'], **kwargs)
        # dataSetXS = dataSetMatE.copy()   # Just set as a copy (redundant here) - gives issues with XS looping below.
        # dataSetXS = None
        dataSetXS = []   # This will skip further processing loop below.

    # Default ePS case
    else:
        # Read fileset
        # 13/10/20 with updated sorting code, this should return
        # - a one-element list for a dir with Eke split files.
        # - a multi-element list for a dir with multiple jobs.
        # - Note cross-over with multiJob class in latter case.
        dataSetXS = readMatEle(fileBase = dataPath, fileIn = fileIn, recordType = 'CrossSection', verbose = self.verbose['sub'])  # Set for XS + betas only
        dataSetMatE = readMatEle(fileBase = dataPath, fileIn = fileIn, recordType = 'DumpIdy', verbose = self.verbose['sub'])

    # Log some details - currently not passed directly from readMatEle()
    # NOTE - with updated code, this is now in data.fileList
#         fList = [item.attrs['fileList'] for item in dataSetXS]
#         fN = len(fList)
#         jobN = len(dataSetXS)

    # Now handled in readMatEle()
    # if len(dataSetXS) > 1:
    #     if self.verbose > 2:
    #         print('Processing data as subdir jobsets, with Eke stacking.')
    #
    #     # Stack multi-E Xarrays into single array.
    #     # Keep results as list for compatibility with rest of code (otherwise will slice Xarray)
    #     dataXS = [xr.combine_nested(dataSetXS, concat_dim = ['Eke']).sortby('Eke')]
    #     dataMatE = [xr.combine_nested(dataSetMatE, concat_dim = ['Eke']).sortby('Eke')]

    # Set dictionary to hold sorted xr.datasets - NOW set in __init__
    dsSet = {}
    jobNotes = []
    fNTotal = 0  # Count total files

    # Set other attribs from source files
    for m, item in enumerate(dataSetXS):
        # Set job info from first datafile of each set (orbital)
        # For Eke stated cases, this will use single file only, and assume identical props.
        dataFile = Path(item.attrs['fileBase'], item.attrs['file'])
        dataSetXS[m].attrs['jobInfo'] = headerFileParse(dataFile, verbose = self.verbose['sub'])
        dataSetXS[m].attrs['molInfo'] = molInfoParse(dataFile, verbose = self.verbose['sub'])

        # Set orb info
        dataSetXS[m].attrs['orbX'], dataSetXS[m].attrs['orbInfo'] = getOrbInfo(item.attrs['jobInfo'], item.attrs['molInfo'])

        # 07/10/22 - skip for empty comments case
        if dataSetXS[m].attrs['jobInfo']['comments']:
            # Additional labels, use these in plotting routines later
            # Try/except here to allow for different formats, fallback to full comment line.
            try:
                dataSetXS[m].attrs['jobLabel'] = item.jobInfo['comments'][1].split('(', maxsplit=1)[1].split(')')[0]
                dataSetMatE[m].attrs['jobLabel'] = item.jobInfo['comments'][1].split('(', maxsplit=1)[1].split(')')[0]
            except IndexError:
                dataSetXS[m].attrs['jobLabel'] = item.jobInfo['comments'][1]
                dataSetMatE[m].attrs['jobLabel'] = item.jobInfo['comments'][1]


        # 07/10/22 - skip for missing IPot case
        if 'IPot' in dataSetXS[m].attrs['jobInfo'].keys():
            # Set absolute photon energy
            dataSetXS[m]['Ehv'] = (item['Ehv'] - (float(item.jobInfo['IPot']) + item.orbX['E'].data[0])).round(self.Edp)
            dataSetMatE[m]['Ehv'] = (dataSetMatE[m]['Ehv'] - (float(item.jobInfo['IPot']) + item.orbX['E'].data[0])).round(self.Edp)

        # jobNotes.append({ 'batch': dataXS[m].jobInfo['comments'][0].strip('#').strip(),
        #                 'event': dataXS[m].jobInfo['comments'][1].split(',', maxsplit=1)[1].strip(),
        #                 'orbLabel': dataXS[m].jobInfo['comments'][1].split('(', maxsplit=1)[1].split(')')[0],
        #                 'orbE': dataXS[m].orbX['E'].data[0]
        #                 })

        # Job notes for plot labels etc.
        # This is basically as per summary.jobInfo() routine
        try:
            jobNotes.append({ 'batch': dataSetXS[m].jobInfo['comments'][0].strip('#').strip(),
                            'event': dataSetXS[m].jobInfo['comments'][1].split(',', maxsplit=1)[1].strip(),
                            'orbLabel': dataSetXS[m].jobInfo['comments'][1].split('(', maxsplit=1)[1].split(')')[0],
                            'orbE': dataSetXS[m].orbX['E'].data[0]
                            })
        # Could do with some proper pattern-matching here, for now just fall-back to [-1] element and hope for the best!
        except IndexError:
            jobNotes.append({ 'batch': dataSetXS[m].jobInfo['comments'][0].strip('#').strip(),
                            'event': dataSetXS[m].jobInfo['comments'][1].split(',', maxsplit=1)[-1].strip('#'),
                            'orbLabel': dataSetXS[m].jobInfo['comments'][1].split('(', maxsplit=1)[-1].split(')')[0],
                            'orbE': dataSetXS[m].orbX['E'].data[0]
                            })

        # Job notes for plot labels etc.
        # This is basically as per summary.jobInfo() routine
        # TODO: Should replace/update/consolidate with regex and more to util func.
        if self.verbose['main'] > 1:
            print(f"Batch: {jobNotes[-1]['batch']}")
            print(f"Orbital: {jobNotes[-1]['orbLabel']}")
            print(f"OrbE: {jobNotes[-1]['orbE']}")

        #*** Set outputs to dict.
        # Set key
        if keyType == 'orb':
            key = f"orb{item.orbX['orb'].data[0]}"

            if key in dsSet.keys():
                key = f"orb{item.orbX['orb'].data[0]}-{m}"  # Add m if key already exists to avoid overwrite.

        # Use int key. This might cause issues for multiJob wrapper
        elif keyType == 'int':
            key = m

        # 06/04/21 Crude hack for multiJob case to pass preset key (for dir stacking with no overwrite for bond scan case)
        else:
            key = keyType

        # Set as xr.ds(), staked by dataType, one per job/orbital
        # Note stacking all jobs can be problematic due to alignment of dims - so use one Dataset per job as model for now.
        #             dsSet[f"orb{item.orbX['orb'].data[0]}"] = xr.Dataset({'XS':dataSetXS[m], 'matE':dataSetMatE[m]})

        # Issues with using xr.Dataset for objects with some different dimensions - use dict instead?
        # May want to use an xr.Dataset for computed properties however? In that case dims should match.
        dsSet[key] = {'XS':dataSetXS[m], 'matE':dataSetMatE[m],
                      'jobNotes':jobNotes[-1]}

        # Set additional metadata, was previously in self.job, but set here for shared key
        if isinstance(item.attrs['fileList'], list):
            fList = item.attrs['fileList']
        else:
            fList = [item.attrs['fileList']]

        fN = len(fList)

        # dsSet[key]['job'] = {'dir': dataPath, 'fN': fN, 'files': fList}
        dsSet[key]['job'] = {'dir': item.attrs['fileBase'], 'fN': fN, 'files': fList}  # Set individual fileBase here for cases when dataPath=None
        fNTotal += fN


    # 07/10/22 - case for no XS data (e.g. different file type)
    if not dataSetXS:
        # Set other attribs from source files
        for m, item in enumerate(dataSetMatE):
            # Set job info from first datafile of each set (orbital)
            # For Eke stated cases, this will use single file only, and assume identical props.
            dataFile = Path(item.attrs['fileBase'], item.attrs['file'])
            dataSetMatE[m].attrs['jobInfo'] = headerFileParse(dataFile, verbose = self.verbose['sub'])
            dataSetMatE[m].attrs['molInfo'] = molInfoParse(dataFile, verbose = self.verbose['sub'])


            #*** Set outputs to dict.
            # Set key
            if keyType == 'orb':
                # key = f"orb{item.orbX['orb'].data[0]}"
                #
                # if key in dsSet.keys():
                #     key = f"orb{item.orbX['orb'].data[0]}-{m}"  # Add m if key already exists to avoid overwrite.

                # Orb not set in this case - default to int.
                key = m

            # Use int key. This might cause issues for multiJob wrapper
            elif keyType == 'int':
                key = m

            # 07/10/22 - use name?
            elif keyType == 'name':
                key = item.name

            # 06/04/21 Crude hack for multiJob case to pass preset key (for dir stacking with no overwrite for bond scan case)
            else:
                key = keyType

            # Set as xr.ds(), staked by dataType, one per job/orbital
            # Note stacking all jobs can be problematic due to alignment of dims - so use one Dataset per job as model for now.
            #             dsSet[f"orb{item.orbX['orb'].data[0]}"] = xr.Dataset({'XS':dataSetXS[m], 'matE':dataSetMatE[m]})

            # Issues with using xr.Dataset for objects with some different dimensions - use dict instead?
            # May want to use an xr.Dataset for computed properties however? In that case dims should match.
            dsSet[key] = {'XS':None, 'matE':dataSetMatE[m],
                          'jobNotes':None}

            # Set additional metadata, was previously in self.job, but set here for shared key
            if isinstance(item.attrs['fileList'], list):
                fList = item.attrs['fileList']
            else:
                fList = [item.attrs['fileList']]

            fN = len(fList)

            # dsSet[key]['job'] = {'dir': dataPath, 'fN': fN, 'files': fList}
            dsSet[key]['job'] = {'dir': item.attrs['fileBase'], 'fN': fN, 'files': fList}  # Set individual fileBase here for cases when dataPath=None
            fNTotal += fN




    # Set self
    self.data.update(dsSet)
    self.jobNotes.append(jobNotes)

    # Meta data - NOW MAINLY MOVED TO data['job']
    # Just use this for some aggregate stuff (per dir)
    # NOTE - this is not persistent over multiple IO cycles, may want to modify?
#         self.job['files'] = fList
    self.job['fN'] = fNTotal
#         self.job['jobN'] = jobN
#         self.job['files']

    # Propagate raw data for testing only.
    # Note this is currently only set for single IO cycle, so will overwrite any existing data
    self.dataSetXS = dataSetXS
    self.dataSetMatE = dataSetMatE

    if self.verbose['main']:
        self.jobsSummary()



def jobsSummary(self):
    """
    Print some general info.

    TODO: add some info!
    """

    # Set pprint object
    pp = pprint.PrettyPrinter(indent=4)


#         fN = [len(item['fList']) for item in self.job['files'].values()]
    # print(f"Dir {self.job['fileBase']}, with {self.job['fN']} files.")
#         print(f"Job structure: {self.jobs['jobStructure']}")

#         if self.jobs['jobStructure'] == 'subDirs':
#             print("(Results E stacked to one dataset per dir.)")

    # Print job details
#         [ep.jobSummary(testClass.dataSets[0]['XS'][1].jobInfo);
    for key in self.data:
        print(f"\n*** Job {key} details")
        print(f"Key: {key}")

        if 'job' in self.data[key].keys():
            print(f"Dir {self.data[key]['job']['dir']}, {self.data[key]['job']['fN']} file(s).")
            pp.pprint(self.data[key]['jobNotes'])
    #             print(f"Directory: {self.jobs['files'][key]['dir']}")
    #             print(f"{self.job['files']['fN']} files")

            if self.verbose['main'] > 1:
                print("File list: ", *self.data[key]['job']['files'], sep='\n   ')

        else:
            print(f"No 'job' info set for self.data[{key}].")

#             for m, item in enumerate(self.data[key]['XS']):
            # [print(line.strip('#')) for line in self.dataSets[key]['XS'][m].jobInfo['comments'][0:4]]
            # print(self.dataSets[key]['jobNotes'][m])
#                 pp.pprint(self.data[key]['jobNotes'])



# Mol info
# Add this, assuming that first job is representative (should add logic to check this)
# Based on existing code in ep.util.jobSummary()
def molSummary(self, dataKey = None, tolConv = 1e-2):

    # Check/get key from index - bit ugly, should decide to pass key or index here?
    # In either case, want default to be 1st entry
    if dataKey is None:
        key = list(self.data.keys())[0]
        dataKey = key


    molInfo = self.data[dataKey]['XS'].attrs['molInfo']

    # Plot molecular structure (crudely)
    print("*** Molecular structure")
    molPlot(molInfo)

    # Print orbTable
    print("\n*** Molecular orbital list (from ePS output file)")
    print("EH = Energy (Hartrees), E = Energy (eV), NOrbGrp, OrbGrp, GrpDegen = degeneracies and corresponding orbital numbering by group in ePS, NormInt = single centre expansion convergence (should be ~1.0).")
#     print(molInfo['orbTable'].to_pandas()) # print() here doesn't give nice interactive result in notebook, although returning table does.
    orbPD = molInfo['orbTable'].to_pandas()

    # Add symmery labels
    # 20/05/22 updated for changed return type.
    # orbPD.insert(1, "Sym", molInfo['orbTable'].Sym) # THIS now adds xarray data types
    orbPD.insert(1, "Sym", molInfo['orbTable'].Sym.values)  # OK
    # orbPD.insert(1, "Sym", molInfo['orbTable'].Sym.to_pandas())  # also OK

    # Tidy a bit...
    orbPD.drop(columns=["N","Type"], inplace=True)

    # Check if notebook and output
    if isnotebook():
        display(orbPD) # Use IPython display(), works nicely for notebook output
        # Q: is display() always loaded automatically? I think so.
    else:
        print(orbPD)

    # Check ExpOrb outputs...
    ind = (molInfo['orbTable'][:,8].values < 1-tolConv) + (molInfo['orbTable'][:,8].values > 1+tolConv)
    if ind.any():
        print(f"\n*** Warning: some orbital convergences outside single-center expansion convergence tolerance ({tolConv}):")
        print(molInfo['orbTable'][ind, [0, 8]].values)

#     return orbPD  # Could also set return value here for simple orbPD printing.


def matEtoPD(self, keys = None, xDim = 'Eke', Erange = None, dataType = 'matE', printTable = True, selDims = None, pType = None,
            thres = None, drop = True, fillna = False, squeeze = True, setPD = True):
    """
    Convert Xarray to PD for nice tabular display.

    Basically code as per basicPlotters.lmPlot(), but looped over datasets.

    30/10/20 Added & reworked from multiJob test code.
            Changed output to nest in existing Xarray & allow multiple datatypes

    """

    keys = self._keysCheck(keys)
    # Default to all datasets
    # if keys is None:
    #     keys = self.dataSets.keys()

    pdConv = [] # List for outputs

    for key in keys:
        # Init empty list for daPlotpd data
        # if setPD:
        #     self.data[key]['daPlotpd'] = []

        # for m, item in enumerate(self.data[key]['matE']):

        #*** Set, select & threshold
        da = self.data[key][dataType]

        # Check xDim Eke/Ehv and clip if specified
        # More elegant way to swap on dims?
        if xDim == 'Ehv':
            # Subset before plot to avoid errors on empty array selection!
            da = da.swap_dims({'Eke':'Ehv'})

        daSub = matEleSelector(da, thres=thres, inds = selDims, dims = xDim, drop = drop)

        # NOTE: assumes xDim = Eke or Ehv, otherwise will fail.
        if Erange is not None:
            # Check Etype, assuming only Eke or Ehv options.
            Etype = 'Eke'
            if 'Ehv' in daSub.dims:
                Etype = 'Ehv'

            daSub = daSub.sel(**{Etype:slice(Erange[0], Erange[1])})  # With dict unpacking for var as keyword

        #*** Data conversion if specified
        if pType is not None:
            daSub = plotTypeSelector(daSub, pType = pType, axisUW = xDim)

        #*** Convert to PD
        daPD, daSub = multiDimXrToPD(daSub, colDims = xDim, thres = thres, dropna = drop, fillna = fillna, squeeze = squeeze)

        # pdConv.append(daPD)

        if printTable:
            print(f"\n*** {da.jobLabel}")
            print(f"Matrix element table, threshold={thres}, data type={daSub.dtype}.")
            # Check if notebook and output
            if isnotebook():
                display(daPD) # Use IPython display(), works nicely for notebook output
                # Q: is display() always loaded automatically? I think so.
            else:
                print(daPD)

        # Set Pandas table to dataset if specified.
        if setPD:
            # self.data[key]['daPlotpd'].append(daPD)
            self.data[key][dataType].attrs['pd'] = daPD


    # Return value if not set to dataSets.
    if not setPD:
        # return pdConv
        return daPD
