"""
Core classes for ePSproc data.

13/10/20    v1  Started class development, reverse engineering a little from multiJob.py case.


"""

# Local functions
from epsproc import readMatEle, headerFileParse, molInfoParse, lmPlot, matEleSelector, plotTypeSelector, multiDimXrToPD, mfpad, sphSumPlotX
# from epsproc.classes.base import ePSbase
from epsproc.util.summary import getOrbInfo, molPlot

class ePSbase():
    """
    Base class for ePSproc.

    Define data model for a single ePS job, defined as a specific ionization event (but may involve multiple ePS output files over a range of Ekes and symmetries).

    Define methods as wrappers for existing ePSproc functions on self.data.

    13/10/20    v1, pulling code from ePSmultiJob().

    - Read datasets from a single dir (uses :py:func:`epsproc.readMatEle`).
    - Sort to dictionaries and Xarray datasets (needs some work).
    - Basic selection, plotting and calculation wrappers in development.


    Parameters
    ----------
    fileBase : str or Path object, default = None
        Base directory to scan for ePS files, subdirs will NOT be searched.
        Use ePSmultiJob class for multi-dir scanning case.

    prefix : str, optional, default = None
        Set prefix string for file checks (cf. wfPlot class).
        Only necessary if automated file sorting fails.

    ext : str, optional, default = '.out'
        Set default file extension for dir scanning.
        This should match the file extension for ePolyScat output files.

    Edp : int, optional, default = 2
        Set default dp for Ehv conversion. May want to set this elsewhere instead... maybe just for plotting?
        TODO: also consider axis reindex, lookups and interp functions here - useful for differences between datasets.

    verbose : int, optional, default = 1
        Set verbosity level for printing/error checking.
        Not yet fully implemented, but, generally:
        - 0, no printed output.
        - 1, basic printed info.
        - 2, print all info, including subfunction outputs.


    TODO:

    - verbosity levels, subtract for subfunctions? Or use a dict to handle multiple levels?


    """

    # TODO: set to varg in for jobs dict
    def __init__(self, fileBase = None,
                 prefix = None, ext = '.out', Edp = 1, verbose = 1):

        self.verbose = verbose

        self.Edp = Edp  # Set default dp for Ehv conversion. May want to set this elsewhere, and match to data resolution?
        self.lmPlotOpts = {}  # Set empty dict for plot options.

        # Quick hack to allow for use of Path() setting below.
        # In this case will use jobDirs instead.
        if fileBase is None:
            fileBase = ''

        # Set file properties
        self.job = {'fileBase':Path(fileBase),
                     'ext':ext,
                     'prefix':prefix,
                     }

         # Set master dictionary to hold sorted datasets
         # Set here to allow persistence over multiple scanFiles() runs from child classes
         self.data = {}
         self.jobNotes = []


    def scanFiles(self, dataPath = None):
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

        """

        # Allow override here for subclasses/independent use
        if dataPath is None:
            dataPath = self.job['fileBase']

        # For Xarray Dataset stacking...
        # dsXS = xr.Dataset()  # Set blank dataset
        # dsMatE = xr.Dataset()

        # Read fileset
        # 13/10/20 with updated sorting code, this should return
        # - a one-element list for a dir with Eke split files.
        # - a multi-element list for a dir with multiple jobs.
        # - Note cross-over with multiJob class in latter case.
        dataSetXS = readMatEle(fileBase = dataPath, recordType = 'CrossSection', verbose = self.verbose)  # Set for XS + betas only
        dataSetMatE = readMatEle(fileBase = dataPath, recordType = 'DumpIdy', verbose = self.verbose)

        # Log some details - currently not passed directly from readMatEle()
        fN = len(dataSetXS)
        fList = [item.attrs['file'] for item in dataSetXS]

        # Now handled in readMatEle()
        # if len(dataSetXS) > 1:
        #     if self.verbose > 2:
        #         print('Processing data as subdir jobsets, with Eke stacking.')
        #
        #     # Stack multi-E Xarrays into single array.
        #     # Keep results as list for compatibility with rest of code (otherwise will slice Xarray)
        #     dataXS = [xr.combine_nested(dataSetXS, concat_dim = ['Eke']).sortby('Eke')]
        #     dataMatE = [xr.combine_nested(dataSetMatE, concat_dim = ['Eke']).sortby('Eke')]


        # Set other attribs from source files
        for m, item in enumerate(dataSetXS):
            # Set job info from first datafile of each set (orbital)
            # For Eke stated cases, this will use single file only, and assume identical props.
            dataFile = Path(item.attrs['fileBase'], item.attrs['file'])
            dataSetXS[m].attrs['jobInfo'] = headerFileParse(dataFile, verbose = self.verbose)
            dataSetXS[m].attrs['molInfo'] = molInfoParse(dataFile, verbose = self.verbose)

            # Set orb info
            dataSetXS[m].attrs['orbX'], dataSetXS[m].attrs['orbInfo'] = getOrbInfo(item.attrs['jobInfo'], item.attrs['molInfo'])

            # Additional labels, use these in plotting routines later
            dataSetXS[m].attrs['jobLabel'] = item.jobInfo['comments'][1].split('(', maxsplit=1)[1].split(')')[0]
            dataSetMatE[m].attrs['jobLabel'] = item.jobInfo['comments'][1].split('(', maxsplit=1)[1].split(')')[0]

            # Set absolute photon energy
            dataSetXS[m]['Ehv'] = (item['Ehv'] - (float(item.jobInfo['IPot']) + item.orbX['E'].data[0])).round(self.Edp)
            dataSetMatE[m]['Ehv'] = (item['Ehv'] - (float(item.jobInfo['IPot']) + item.orbX['E'].data[0])).round(self.Edp)

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
            if self.verbose > 1:
                print(f"Batch: {jobNotes[-1]['batch']}")
                print(f"Orbital: {jobNotes[-1]['orbLabel']}")
                print(f"OrbE: {jobNotes[-1]['orbE']}")

            #*** Set outputs to dict.
            # Set key
            key = f"orb{item.orbX['orb'].data[0]}"

            if key in dsSet.keys():
                key = f"orb{item.orbX['orb'].data[0]}-{m}"  # Add m if key already exists to avoid overwrite.

            # Set as xr.ds(), staked by dataType, one per job/orbital
            # Note stacking all jobs can be problematic due to alignment of dims - so use one Dataset per job as model for now.
            #             dsSet[f"orb{item.orbX['orb'].data[0]}"] = xr.Dataset({'XS':dataSetXS[m], 'matE':dataSetMatE[m]})

            # Issues with using xr.Dataset for objects with some different dimensions - use dict instead?
            # May want to use an xr.Dataset for computed properties however? In that case dims should match.
            dsSet[key] = {'XS':dataSetXS[m], 'matE':dataSetMatE[m],
                            'jobNotes':jobNotes[-1]}

        # Set self
        self.data.update(dsSet)
        self.jobNotes.append(jobNotes)

        # Propagate raw data for testing only.
        self.dataSetXS = dataSetXS
        self.dataSetMatE = dataSetMatE

        if self.verbose:
            self.jobsSummary()
