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

    Define data for a single ePS job, defined as a specific ionization event (but may involve multiple ePS output files over a range of Ekes and symmetries).

    Define methods as wrappers for existing ePSproc functions on self.data.

    13/10/20    v1, pulling code from ePSmultiJob().

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

    def scanFiles(self, dataPath = None):

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

        # Set dictionary to hold sorted xr.datasets
        dsSet = {}
        jobNotes = []

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



            # Set as xr.ds(), staked by dataType, one per job/orbital
            # Note stacking all jobs can be problematic due to alignment of dims - so use one Dataset per job as model for now.
            #             dsSet[f"orb{item.orbX['orb'].data[0]}"] = xr.Dataset({'XS':dataSetXS[m], 'matE':dataSetMatE[m]})

            # Issues with using xr.Dataset for objects with some different dimensions - use dict instead?
            # May want to use an xr.Dataset for computed properties however? In that case dims should match.
            dsSet[f"orb{item.orbX['orb'].data[0]}"] = {'XS':dataSetXS[m], 'matE':dataSetMatE[m], 'jobNotes':jobNotes[-1]}

        # Set self
        self.data = dsSet
        self.jobNotes = jobNotes

        # Propagate raw data for testing only.
        self.dataSetXS = dataSetXS
        self.dataSetMatE = dataSetMatE
