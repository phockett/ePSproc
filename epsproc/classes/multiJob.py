"""
Core classes for ePSproc data to handle multiple job filesets (energy and/or orbitals etc.)

13/10/20    V2  Consolidating with ePSbase class, some functionality moved there.

23/09/20    Added ePSmultiJob class, currently very rough.

14/09/20    v1  Started class development. See ePSproc_multijob_class_dev_140920_bemo.ipynb

"""

from pathlib import Path
from matplotlib import pyplot as plt  # For plot legends with Xarray plotter

import numpy as np # Needed only for np.nan at the moment
import xarray as xr
import pprint

# Local functions
from epsproc import readMatEle, headerFileParse, molInfoParse, lmPlot, matEleSelector, plotTypeSelector, multiDimXrToPD, mfpad, sphSumPlotX
from epsproc.classes.base import ePSbase
from epsproc.util.summary import getOrbInfo, molPlot
from epsproc.util.env import isnotebook
from epsproc.plot import hvPlotters

# Class for multiple ePS datasets
class ePSmultiJob(ePSbase):
    """
    Class for handling multiple ePS jobs/datasets.

    - Read datasets from a single dir, or set of dirs.
    - Sort to dictionaries and Xarray datasets (needs some work).
    - Comparitive plotting (in development).


    Parameters
    ----------
    fileBase : str or Path object, default = None
        Base directory to scan for ePS files, subdirs will also be searched.
        This is required unless jobDirs is set instead.

    jobDirs : list of str or Path objects, default= None
        List of dirs containing ePS output files to read, subdirs will NOT be searched.

    jobStructure : str, optional, default = None
        This will be set automatically by self.scanDirs(), but can also be passed to override.
        Values "subDirs" or "rootDir", in former case multiple files will be stacked by Eke.

    prefix : str, optional, default = None
        Set prefix string for file checks (cf. wfPlot class). NOT YET USED.

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
        UPDATE: now handled in base class with verbose['main'] and verbose['sub'].



    """

    # TODO: set to varg in for jobs dict
    def __init__(self, fileBase = None, jobDirs = None, jobStructure = None,
                 prefix = None, ext = '.out', Edp = 1, verbose = 1):

        # self.verbose = verbose
        #
        # self.Edp = Edp  # Set default dp for Ehv conversion. May want to set this elsewhere, and match to data resolution?
        # self.lmPlotOpts = {}  # Set empty dict for plot options.

        # Quick hack to allow for use of Path() setting below.
        # In this case will use jobDirs instead.
        if fileBase is None:
            fileBase = ''

        # Run rest of __init__ from base class
        super().__init__(fileBase = fileBase, prefix = prefix, ext = ext, Edp = Edp, verbose = verbose)

        # Set additional properties for multiJob case
        # Set file properties
        self.jobs = {'fileBase':Path(fileBase),
                     'ext':ext,
                     'prefix':prefix,
                     'jobDirs':jobDirs,
                     'jobStructure':jobStructure}

        if jobDirs is None:
            if verbose:
                print(f"Scanning {fileBase} for ePS jobs.")

            self.scanDirs()


    def scanDirs(self):
        """
        Scan dir structure for folders containing ePS files.

        Compatibility... this assumes one of two dir structures:

        - Old structure, with multiple job output files per dir (but not per E range).
        - New structure, with multiple output files per E range, subdirs per job/orb.

        This is set in jobStructure variable, as "rootDir" or "subDirs" respectively.

        """

        # Scan dirs with pathlib
        # https://stackoverflow.com/a/30925692, + glob for .out files
        # Subdirs only - NOT ROOT!
        testSub = [p for p in self.jobs['fileBase'].iterdir() if (p.is_dir() and bool(list(p.glob(f"*{self.jobs['ext']}"))))]

        if testSub and (self.jobs['jobStructure'] is None or 'subDirs'):
            self.jobs['jobDirs'] = testSub
            self.jobs['jobStructure'] = 'subDirs'

            if self.verbose['main'] > 0:
                print(f"Found ePS output files in subdirs: {testSub}")


        # Check root dir
        testRoot = list(self.jobs['fileBase'].glob(f"*{self.jobs['ext']}"))

        if testRoot and (self.jobs['jobStructure'] is None or 'rootDir'):
            self.jobs['jobDirs'] = [self.jobs['fileBase']]
            self.jobs['jobStructure'] = 'rootDir'

            if self.verbose['main'] > 0:
                print(f"Found ePS output files in root dir {self.jobs['fileBase']}.")


    def scanFiles(self, keys = None, outputKeyType = 'dir'):
        """
        Scan ePS output files from dir list.

        Adapted from https://phockett.github.io/ePSdata/XeF2-preliminary/XeF2_multi-orb_comparisons_270320-dist.html

        Currently outputting complicated dataSets dictionary/list - need to sort this out!
        - Entry per dir scanned.
        - Sub-entries per file, but collapsed in some cases.
        - Sending to Xarray dataset with orb labels should be cleaner.

        Parameters
        -----------
        keys : int or list, optional, defaults to all jobDirs
            Set keys (dirs) to use by index.
            Default = enumerate(self.jobs['jobDirs'])


        outputKeyType : str, optional, default = 'orb'
            Types as supported by super().scanFiles()
            'orb': Use orbital labels as dataset keys
            'int': Use integer labels as dataset keys (will be ordered by file read)
            Any other setting will result in key = keyType, which can be used to explicitly pass a key (e.g. in multijob wrapper case). This should be tidied up.



        TODO:

        - Flatten data structure to remove unnecessary nesting.
        - Fix keys to apply to single dir case (above should also fix this).
        - convert outputs to Xarray dataset. Did this before, but currently missing file (on AntonJr)! CHECK BACKUPS - NOPE.
        - Confirm HV scaling - may be better to redo this, rather than correct existing values?

        - Fix xr.dataset: currently aligns data, so will set much to Nan if, e.g., different symmetries etc.
        Change to structure as ds('XS','matE') per orb, rather than ds('XS') and ds('matE') for all orbs?
        This should also be in line with hypothetical base dataclass, which will be per orb by defn.
        UPDATE 16/10/20: now handled in base class, use flat dict with entry per job or dir (for E-stacked case).


        14/10/20    v2  Now using ePSbase() class for functionality.
                        DATASTRUCTURE IS NOW FLAT

        """

        # Original version
        dataSets = {}  # Set dict to hold all data for the moment.
        jobs = {}  # Set dict to hold summary info

        # For Xarray Dataset stacking...
        dsXS = xr.Dataset()  # Set blank dataset
        dsMatE = xr.Dataset()

        # 29/09/20 Quick hack to allow subselection by index.
        # Should rethink this and add some better selection options (see wf IO code...?)
        if keys is None:
            keys = enumerate(self.jobs['jobDirs'])
        else:
            if isinstance(keys, int):
                keys = [keys]

            keys = [(n, self.jobs['jobDirs'][n]) for n in keys]

        for n, dirScan in keys:
#             dataPath = Path(workingDir[0], dirScan)
            dataPath = dirScan  # Assume abs paths (?)

            # 06/04/21 Crude hack for multiJob case to pass preset key (for dir stacking with no overwrite for bond scan case)
            if outputKeyType == 'dir':
                outputKey = dirScan.name
            else:
                outputKey = outputKeyType

            # Call super class to handle reading per dir
            super().scanFiles(dataPath = dataPath, reset = False, keyType = outputKey)
            #********** IN PROGRESS!

#             # For dir scan
#             dataSetXS = readMatEle(fileBase = dataPath, recordType = 'CrossSection', verbose = self.verbose)  # Set for XS + betas only
#             dataSetMatE = readMatEle(fileBase = dataPath, recordType = 'DumpIdy', verbose = self.verbose)
#
#             # Log some details - currently not passed directly from readMatEle()
#             fN = len(dataSetXS)
            # fList = [item.attrs['file'] for item in dataSetXS]

            # fN = len(self.dataSetXS)
            # fList = [item.attrs['file'] for item in self.dataSetXS]

#             if self.jobs['jobStructure'] == 'subDirs' and len(dataSetXS) > 1:
#                 if self.verbose['main'] > 2:
#                     print('Processing data as subdir jobsets, with Eke stacking.')
#
#                 # Stack multi-E Xarrays into single array.
#                 # Keep results as list for compatibility with rest of code (otherwise will slice Xarray)
#                 dataXS = [xr.combine_nested(dataSetXS, concat_dim = ['Eke']).sortby('Eke')]
#                 dataMatE = [xr.combine_nested(dataSetMatE, concat_dim = ['Eke']).sortby('Eke')]
#
#                 # Set job info from first datafile of each set (orbital)
#                 dataFile = Path(dataXS[0].attrs['fileBase'], dataXS[0].attrs['file'])
#                 dataXS[0].attrs['jobInfo'] = headerFileParse(dataFile, verbose = self.verbose)
#                 dataXS[0].attrs['molInfo'] = molInfoParse(dataFile, verbose = self.verbose)
#
#                 # Set orb info
#                 dataXS[0].attrs['orbX'], dataXS[0].attrs['orbInfo'] = getOrbInfo(dataXS[0].attrs['jobInfo'], dataXS[0].attrs['molInfo'])
#
#                 # Additional labels, use these in plotting routines later
#                 dataXS[0].attrs['jobLabel'] = dataXS[0].jobInfo['comments'][1].split('(', maxsplit=1)[1].split(')')[0]
#                 dataMatE[0].attrs['jobLabel'] = dataXS[0].jobInfo['comments'][1].split('(', maxsplit=1)[1].split(')')[0]
#
#                 # Set absolute photon energy
#                 dataXS[0]['Ehv'] = (dataXS[0]['Ehv'] - (float(dataXS[0].jobInfo['IPot']) + dataXS[0].orbX['E'].data[0])).round(self.Edp)
#                 dataMatE[0]['Ehv'] = (dataMatE[0]['Ehv'] - (float(dataXS[0].jobInfo['IPot']) + dataXS[0].orbX['E'].data[0])).round(self.Edp)
#
#                 # Set as xr.ds(), label by orbital
#                 # TODO: fix this, in some cases gives errors with multiple values - probably an Ehv issue?
#                 # dsXS[dataXS[0].orbX['orb'].data[0]] = dataXS[0]
#                 # dsMatE[dataXS[0].orbX['orb'].data[0]] = dataMatE[0]
#                 # With str labelling, as per non-subdir case - was this the issue?
#                 # dsXS[f"orb{dataXS[0].orbX['orb'].data[0]}"] = dataXS[0]
#                 # dsMatE[f"orb{dataXS[0].orbX['orb'].data[0]}"] = dataMatE[0]
#                 # Main issue was aligning to coords, e.g. need xr.concat([data1, data2], "Sym").to_dataset() in general to properly stack.
#
#                 # Job notes for plot labels etc.
#                 # This is basically as per summary.jobInfo() routine
#                 # TODO: Should replace/update/consolidate with regex and more to util func.
#                 if self.verbose['main'] > 1:
#                     print(f"Batch: {dataXS[0].jobInfo['comments'][0].strip('#').strip()}")
#                     print(f"Orbital: {dataXS[0].jobInfo['comments'][1].split(',', maxsplit=1)[1].strip()}")
#                     print(f"OrbE: {dataXS[0].orbX['E'].data[0]}")
#
#                 jobNotes = []
#                 jobNotes.append({ 'batch': dataXS[0].jobInfo['comments'][0].strip('#').strip(),
#                                 'event': dataXS[0].jobInfo['comments'][1].split(',', maxsplit=1)[1].strip(),
#                                 'orbLabel': dataXS[0].jobInfo['comments'][1].split('(', maxsplit=1)[1].split(')')[0],
#                                 'orbE': dataXS[0].orbX['E'].data[0]
#                                 })
#
#             else:
#                 if self.verbose['main'] > 2:
#                     print('Processing data as single dir jobset, no Eke stacking.')
#
#                 dataXS = dataSetXS
#                 dataMatE = dataSetMatE
#                 jobNotes = []
#
#                 # Set job info for each file
#                 for m, item in enumerate(dataXS):
#
#                     dataFile = Path(item.attrs['fileBase'], dataXS[m].attrs['file'])
#                     dataXS[m].attrs['jobInfo'] = headerFileParse(dataFile, verbose = self.verbose)
#                     dataXS[m].attrs['molInfo'] = molInfoParse(dataFile, verbose = self.verbose)
#
#                     # Set orb info
#                     dataXS[m].attrs['orbX'], dataXS[m].attrs['orbInfo'] = getOrbInfo(dataXS[m].attrs['jobInfo'], dataXS[m].attrs['molInfo'])
#
#                     # Additional labels, use these in plotting routines later
#                     dataXS[m].attrs['jobLabel'] = dataXS[m].jobInfo['comments'][1].split('(', maxsplit=1)[1].split(')')[0]
#                     dataMatE[m].attrs['jobLabel'] = dataXS[m].jobInfo['comments'][1].split('(', maxsplit=1)[1].split(')')[0]
#
#                     # Set absolute photon energy
#                     dataXS[m]['Ehv'] = (dataXS[m]['Ehv'] - (float(dataXS[m].jobInfo['IPot']) + dataXS[m].orbX['E'].data[0])).round(self.Edp)
#                     dataMatE[m]['Ehv'] = (dataMatE[m]['Ehv'] - (float(dataXS[m].jobInfo['IPot']) + dataXS[m].orbX['E'].data[0])).round(self.Edp)
#
#                     # Set as xr.ds(), label by orbital
#                     dsXS[f"orb{dataXS[m].orbX['orb'].data[0]}"] = dataXS[m]
#                     dsMatE[f"orb{dataXS[m].orbX['orb'].data[0]}"] = dataMatE[m]
#
#
#                     # Job notes for plot labels etc.
#                     # This is basically as per summary.jobInfo() routine
#                     try:
#                         jobNotes.append({ 'batch': dataXS[m].jobInfo['comments'][0].strip('#').strip(),
#                                         'event': dataXS[m].jobInfo['comments'][1].split(',', maxsplit=1)[1].strip(),
#                                         'orbLabel': dataXS[m].jobInfo['comments'][1].split('(', maxsplit=1)[1].split(')')[0],
#                                         'orbE': dataXS[m].orbX['E'].data[0]
#                                         })
#                     # Could do with some proper pattern-matching here, for now just fall-back to [-1] element and hope for the best!
#                     except IndexError:
#                         jobNotes.append({ 'batch': dataXS[m].jobInfo['comments'][0].strip('#').strip(),
#                                         'event': dataXS[m].jobInfo['comments'][1].split(',', maxsplit=1)[-1].strip('#'),
#                                         'orbLabel': dataXS[m].jobInfo['comments'][1].split('(', maxsplit=1)[-1].split(')')[0],
#                                         'orbE': dataXS[m].orbX['E'].data[0]
#                                         })
#
# #                     print(m)
# #                     print(dataXS[m].orbX['orb'].data[0])
#
# #             # Set job info from first datafile of each set (orbital)
# #             dataFile = Path(dataXS[0].attrs['fileBase'], dataXS[0].attrs['file'])
# #             dataXS[0].attrs['jobInfo'] = ep.headerFileParse(dataFile)
# #             dataXS[0].attrs['molInfo'] = ep.molInfoParse(dataFile)
# #             dataXS.attrs['jobInfo'] = ep.headerFileParse(dataFile)
# #             dataXS.attrs['molInfo'] = ep.molInfoParse(dataFile)
#
#             dataSets[n] = {}
#             dataSets[n]['dir'] = dirScan
#             dataSets[n]['XS'] = dataXS
#             dataSets[n]['matE'] = dataMatE
#             dataSets[n]['jobNotes'] = jobNotes

            # Dir aggregate data
            # Set from self.job, this data is overwritten for each dir
            # jobs[n] = {'dir': dirScan, 'fN': self.job['fN'], 'fList': self.job['files']}
            jobs[n] = {'dir': dirScan, 'fN':self.job['fN'],
                        'fileList': [self.data[key]['job']['files'] for key in self.data.keys()]}

        self.dataSets = self.data  # Just set pointer for now, but should be redundant in reworked class with flat data dict

        # Aggregate data per dir, from self.data[key][job]
        self.jobs['files'] = jobs  # self.dataSets['jobs']

        # xr.ds types
        # self.dsXS = dsXS
        # self.dsMatE = dsMatE

        if self.verbose['main']:
            self.jobsSummary()

    def jobLabel(self, key = None, lString = None, append=True):
        """
        Reset or append text to jobLabel.
        Very basic.

        TODO: consistency over [m], jobLabel vs. orbLabel rationalisation.

        """

        if append:
            lString = self.dataSets[key]['XS'][0].attrs['jobLabel'] + " " + lString

        self.dataSets[key]['XS'][0].attrs['jobLabel'] = lString
        self.dataSets[key]['matE'][0].attrs['jobLabel'] = lString
        self.dataSets[key]['jobNotes']['orbLabel'] = lString


    def jobsSummary(self):
        """
        Print some general info.

        TODO: add some info!
        """

        # Print multi-dir aggregate data
        # [len(item['fList']) for item in self.jobs['files'].values()]
        try:
            fN = [self.jobs['files'][key]['fN'] for key in self.jobs['files'].keys()]
            print(f"Found {len(self.jobs['files'])} directories, with {sum(fN)} files.")

        except KeyError as E:
            print(f"Skipping jobs summary, missing key: {E}. Try running self.jobsSummary() after file IO completes.")

        # Run parent fn - this is OK, but only outputs self.job info at the moment
        super().jobsSummary()

#         # Set pprint object
#         pp = pprint.PrettyPrinter(indent=4)
#
#
#         fN = [len(item['fList']) for item in self.jobs['files'].values()]
#         print(f"Found {len(self.jobs['files'])} directories, with {sum(fN)} files.")
#         print(f"Job structure: {self.jobs['jobStructure']}")
#
#         if self.jobs['jobStructure'] == 'subDirs':
#             print("(Results E stacked to one dataset per dir.)")
#
#         # Print job details
# #         [ep.jobSummary(testClass.dataSets[0]['XS'][1].jobInfo);
#         for key in self.dataSets:
#             print(f"\n*** Job dir {key} details")
#             print(f"Directory: {self.jobs['files'][key]['dir']}")
#             print(f"{self.jobs['files'][key]['fN']} files")
#
#             for m, item in enumerate(self.dataSets[key]['XS']):
#                 # [print(line.strip('#')) for line in self.dataSets[key]['XS'][m].jobInfo['comments'][0:4]]
#                 # print(self.dataSets[key]['jobNotes'][m])
#                 pp.pprint(self.dataSets[key]['jobNotes'][m])




    # Define photon energy scale - now set at read-in.
#     def setHV(self):

#         # Set Ehv scale - currently incorrectly set by 1st vertical IP
#         for n, key in enumerate(self.dataSets):
#         #     dataSets[key][0]['Ehv'] -= (12.35 + orbIPs[n])
#             dataSets[key]['XS'][0]['Ehv'] = dataSets[key]['XS'][0]['Ehv'] - (12.35 + orbIPs[n])
#             dataSets[key]['matE'][0]['Ehv'] = dataSets[key]['matE'][0]['Ehv'] - (12.35 + orbIPs[n])

#             jobInfo['IPot']


# 14/10/20 Moved to parent fn.
#     def plotGetCro(self, pType = 'SIGMA', Erange = None, Etype = 'Eke', keys = None, backend = 'mpl'):
#         """
#         Basic GetCro (cross-section) data plotting for multijob class. Run self.plot.line(x=Etype, col='Type') for each dataset.
#         (See :py:func:`epsproc.classes.ePSmultiJob.plotGetCroComp()` for comparitive plots over datasets.)
#
#         Note this is for LF averaged parameters, for more details see the `ePS starter notes <https://epsproc.readthedocs.io/en/latest/ePS_ePSproc_tutorial/ePS_tutorial_080520.html#Results>`_ for more details.
#
#         Parameters
#         ----------
#         pType : str, optional, default = 'SIGMA'
#             Set data for plotting, either 'SIGMA' (cross-section) or 'BETA' (B2 parameter).
#             If backend = 'hv' this parameter is not used.
#
#         Erange : list of int or float, optional, default = None
#             Set plot range [Emin, Emax]. Defaults to full data range if not set
#
#         Etype : str, optional, default = 'Eke'
#             Set plot dimension, either 'Eke' (electron kinetic energy) or 'Ehv' (photon energy).
#
#         keys : list, optional, default = None
#             Keys for datasets to plot.
#             If None, all datasets will be plotted.
#
#         backend : str, optional, default = 'mpl'
#             Set plotter to use.
#
#             - 'mpl' : Use Matplotlib/native Xarray plotter
#             - 'hv' : use Holoviews via :py:func:`epsproc.plotters.hvPlotters.XCplot()`
#
#         """
#         # Default to all datasets
#         if keys is None:
#             keys = self.dataSets.keys()
#
# #         if self.jobs['jobStructure'] == 'subDirs':
#         for key in keys:
# #                 testClass.dataSets[key]['XS'][0].sel(XC='SIGMA', Eke=slice(Erange[0], Erange[1])).plot.line(x='Eke', col='Type')   # This works
#             for m, item in enumerate(self.dataSets[key]['XS']):
#
#                 # Set default to full range, same for all cases
#                 if Erange is None:
#                     Erange = [self.dataSets[key]['XS'][m][Etype].min().data, self.dataSets[key]['XS'][m][Etype].max().data]
#
#                 subset = self.dataSets[key]['XS'][m]
#                 jobLabel = self.dataSets[key]['XS'][m].jobLabel
#
#                 # More elegant way to swap on dims?
#                 if Etype == 'Ehv':
#                 # Subset before plot to avoid errors on empty array selection!
#                     subset = subset.swap_dims({'Eke':'Ehv'}) #.sel(**{Etype:slice(Erange[0], Erange[1])})   # With dict unpacking for var as keyword
#
#                     # if subset.any():
#                     #     subset.plot.line(x=Etype, col='Type')
#
#                 subset = subset.sel(**{Etype:slice(Erange[0], Erange[1])})   # With dict unpacking for var as keyword
#
#                 if subset.any():
#                     # Propagate attrs for plot labelling
#                     subset.attrs = self.dataSets[key]['XS'][m].attrs
#
#                     if backend == 'mpl':
#                         subset.sel(XC=pType).plot.line(x=Etype, col='Type')
#
#                         if pType == 'SIGMA':
#                             plt.ylabel('XS/Mb')
#                         else:
#                             plt.ylabel(r"$\beta_{LM}$")
#
#                     if backend == 'hv':
#                         layout, *_ = hvPlotters.XCplot(subset, kdims = Etype)
#
#                         print(f"\n*** {jobLabel}")
#                         # Check if notebook and output
#                         if isnotebook():
#                             display(layout) # Use IPython display(), works nicely for notebook output
#                             # Q: is display() always loaded automatically? I think so.
#                         else:
#                             layout  # Not yet tested - not sure what the options are here for hv/Bokeh.
#                                     # May just want to return this object (hv layout)?
#
#                         # Setting title here doesn't work, now set in XCplot()
#                         # display(layout.opts(title=jobLabel))
#
# #         else:

#     # Basically as per plotGetCro, but subselect and put on single plot.
#     # Should do all of this with Holoviews...!
#     def plotGetCroComp(self, pType='SIGMA', pGauge='L', pSym=('All','All'), Erange = None, Etype = 'Eke',  keys = None):
#         """
#         Basic GetCro (cross-section) data plotting for multijob class, comparitive plots.
#         Run self.plot.line(x=Etype) for each dataset after subselection on Gauge and Symmetry, and use single axis.
#         (See :py:func:`epsproc.classes.ePSmultiJob.plotGetCro()` for plots per dataset.)
#
#         Note this is for LF averaged parameters, for more details see the `ePS starter notes <https://epsproc.readthedocs.io/en/latest/ePS_ePSproc_tutorial/ePS_tutorial_080520.html#Results>`_ for more details.
#
#         Parameters
#         ----------
#         pType : str, optional, default = 'SIGMA'
#             Set data for plotting, either 'SIGMA' (cross-section) or 'BETA' (B2 parameter).
#
#         pGauge : str, optional, default = 'L'
#             Set gauge, either 'L' (Length), 'V' (Velocity) or 'M' (Mixed)
#
#         pSym : tuple of strs, optional, default = ('All','All')
#             Select symmetry, (Cont, Targ).
#             Default value will plot all allowed symmetries.
#
#         Erange : list of int or float, optional, default = None
#             Set plot range [Emin, Emax]. Defaults to full data range if not set
#
#         Etype : str, optional, default = 'Eke'
#             Set plot dimension, either 'Eke' (electron kinetic energy) or 'Ehv' (photon energy).
#
#         keys : list, optional, default = None
#             Keys for datasets to plot.
#             If None, all datasets will be plotted.
#
#         """
#
#
#
#         # Comparison plots over orbs
#
#         # from matplotlib import pyplot as plt  # For legend
#         lText = []
#
#         # Default to all datasets
#         if keys is None:
#             keys = self.dataSets.keys()
#
#         for key in keys:
# #                 testClass.dataSets[key]['XS'][0].sel(XC='SIGMA', Eke=slice(Erange[0], Erange[1])).plot.line(x='Eke', col='Type')   # This works
#             for m, item in enumerate(self.dataSets[key]['XS']):
#
#                 # Set default to full range, same for all cases
#                 if Erange is None:
#                     Erange = [self.dataSets[key]['XS'][m][Etype].min().data, self.dataSets[key]['XS'][m][Etype].max().data]
#
#
#                 # More elegant way to swap on dims?
#                 if Etype == 'Ehv':
#                     # Subset before plot to avoid errors on empty array selection!
#                     subset = self.dataSets[key]['XS'][m].swap_dims({'Eke':'Ehv'}).sel(XC=pType, Type=pGauge, Sym=pSym, **{Etype:slice(Erange[0], Erange[1])})  # With dict unpacking for var as keyword
#
#                     # if subset.any():
#                     #     pltObj = subset.plot.line(x=Etype)
#
#                 else:
#                     subset = self.dataSets[key]['XS'][m].sel(XC=pType, Type=pGauge, Sym=pSym, **{Etype:slice(Erange[0], Erange[1])})   # With dict unpacking for var as keyword
#
#                 if subset.any():
#                     pltObj = subset.plot.line(x=Etype)
#                     lText.append(self.dataSets[key]['jobNotes'][m]['orbLabel'])
#
#                 # Label with orb_sym
# #                 lText.append(self.dataSets[key]['XS'][m].attrs['fileBase'].rsplit('/',maxsplit=1)[0])
#
#                 # lText.append(f"Orb {self.dataSets[key]['XS'][m].attrs['orbInfo']['orbN']} ({self.dataSets[key]['XS'][m].attrs['orbInfo']['orbSym'][0]})")
#
#                 # lText.append(self.dataSets[key]['jobNotes'][m]['orbLabel'])
#
#         plt.legend(lText)
#
#         if pType == 'SIGMA':
#             plt.ylabel('XS/Mb')
#         else:
#             plt.ylabel(r"$\beta_{LM}$")

# 14/10/20 Moved to base class
    # def lmPlot(self, Erange = None, Etype = 'Eke', keys = None, refDataKey = None, reindexTol = 0.5, reindexFill = np.nan, setPD = True, **kwargs):
    #     """
    #     Wrapper for :py:func:`epsproc.lmPlot` for multijob class. Run lmPlot() for each dataset.
    #
    #     Parameters
    #     ----------
    #     Erange : list of int or float, optional, default = None
    #         Set plot range [Emin, Emax]. Defaults to full data range if not set
    #
    #     Etype : str, optional, default = 'Eke'
    #         Set plot dimension, either 'Eke' (electron kinetic energy) or 'Ehv' (photon energy).
    #
    #     keys : list, optional, default = None
    #         Keys for datasets to plot.
    #         If None, all datasets will be plotted.
    #
    #     refDataKey : tuple (key,m), optional, default = None
    #         If set, calculate difference plots against reference dataset.
    #         TODO: implement difference plots.
    #         TODO: implement testing logic, may fail without E-axis forcing, and sym summation?
    #
    #     reindexTol : float, optional, default = 0.1
    #         If computing difference data, the reference data is reindexed to ensure E grid matching.
    #         This specifies tolerance (in E units, usually eV) for reindexing.
    #         If this fails, difference plot may be null.
    #
    #     reindexFill : int or float, optional, default = NaN
    #         Value to use for missing values upon reindexing.
    #         Default matches [Xarray.reindex default](http://xarray.pydata.org/en/stable/generated/xarray.DataArray.reindex.html), i.e. NaN, but this may give issues in some cases.
    #
    #     setPD : bool, optional, default = True
    #         Set Pandas array in main dataset?
    #
    #     kwargs : dict, optional, default = {}
    #         Plotting options to pass to :py:func:`epsproc.lmPlot`.
    #         These will also be set in self.lmPlotOpts for further use.
    #         Note that any existing options in self.lmPlotOpts will also be used, or overwritten if matching keys are found.
    #
    #
    #     Notes
    #     -----
    #     Basic scheme from ePSmultijob.plotGetCro, which loops and switches on Eke/Ehv. Should tidy up at some point.
    #
    #     """
    #
    #     # Default to all datasets
    #     if keys is None:
    #         keys = self.dataSets.keys()
    #
    #     # Set lmPlotOpts
    #     # Check passed args vs. self.lmPlotOpts and overwrite
    #     if kwargs:
    #         for key, value in kwargs.items():
    #             self.lmPlotOpts[key] = value
    #
    #     # Set default to full range of 1st dataset, keep same for all cases
    #     # TODO: set per dataset?
    #     if Erange is None:
    #         Erange = [self.dataSets[keys[0]]['matE'][0][Etype].min().data, self.dataSets[keys[0]]['matE'][0][Etype].max().data]
    #
    #     # Set ref dataset if required
    #     if refDataKey is not None:
    #         refData = self.dataSets[refDataKey[0]]['matE'][refDataKey[1]]
    #
    #         if Etype == 'Ehv':
    #             refData = refData.swap_dims({'Eke':'Ehv'})
    #
    #         refData = refData.sel(**{Etype:slice(Erange[0], Erange[1])})
    #         refData.attrs = self.dataSets[refDataKey[0]]['matE'][refDataKey[1]].attrs # Propagate atrrs.
    #
    #     else:
    #         refData = None
    #
    #
    #     # Loop over datasets
    #     for key in keys:
    # #                 testClass.dataSets[key]['XS'][0].sel(XC='SIGMA', Eke=slice(Erange[0], Erange[1])).plot.line(x='Eke', col='Type')   # This works
    #
    #         # Init empty list for daPlotpd data
    #         if setPD:
    #             self.dataSets[key]['daPlotpd'] = []
    #
    #         for m, item in enumerate(self.dataSets[key]['matE']):
    #
    #             # More elegant way to swap on dims?
    #             if Etype == 'Ehv':
    #             # Subset before plot to avoid errors on empty array selection!
    #                 subset = self.dataSets[key]['matE'][m].swap_dims({'Eke':'Ehv'}).sel(**{Etype:slice(Erange[0], Erange[1])})   # With dict unpacking for var as keyword
    #
    #             else:
    #                 subset = self.dataSets[key]['matE'][m].sel(**{Etype:slice(Erange[0], Erange[1])})   # With dict unpacking for var as keyword
    #
    #             # Difference data
    #             # NOTE: ref. data is reindexed here to ensure E-point subtraction (otherwise mismatch will give null results)
    #             # TODO: may want some selection logic here, this may produce empty results in many cases.
    #             if refData is not None:
    #                 subset = subset - refData.reindex(**{Etype:subset[Etype]}, method='nearest', tolerance=reindexTol, fill_value=reindexFill)
    #                 subset.attrs = self.dataSets[key]['matE'][m].attrs  # Propagate attribs
    #                 subset.attrs['jobLabel'] = f"{subset.attrs['jobLabel']} diff with {refData.attrs['jobLabel']}"
    #
    #                 # Slightly crude check for empty result.
    #                 # Note this is different from subset.any() below, which only checks for 0
    #                 # See https://numpy.org/doc/stable/reference/generated/numpy.any.html
    #                 if subset.max().isnull():
    #                     print("*** Warning, difference array is Null. Try a larger reindexTol value.")
    #
    #             # Run lmPlot
    #             if subset.any():
    #                 daPlot, daPlotpd, legendList, gFig = lmPlot(subset, xDim = Etype, **self.lmPlotOpts)
    #             else:
    #                 daPlotpd = None
    #
    #             # Set Pandas table to dataset if specified.
    #             if setPD:
    #                 self.dataSets[key]['daPlotpd'].append(daPlotpd)  # Set to include None cases to keep indexing. Should set as dict instead?


# 14/10/20 now set in base class.
    # # Mol info
    # # Add this, assuming that first job is representative (should add logic to check this)
    # # Based on existing code in ep.util.jobSummary()
    # def molSummary(self, dataKey = None, tolConv = 1e-2):
    #
    #     # Check/get key from index - bit ugly, should decide to pass key or index here?
    #     # In either case, want default to be 1st entry
    #     if dataKey is None:
    #         key = list(data.dataSets.keys())[0]
    #         dataKey = [key, 0]
    #
    #
    #     molInfo = self.dataSets[dataKey[0]]['XS'][dataKey[1]].attrs['molInfo']
    #
    #     # Plot molecular structure (crudely)
    #     print("*** Molecular structure")
    #     molPlot(molInfo)
    #
    #     # Print orbTable
    #     print("\n*** Molecular orbital list (from ePS output file)")
    #     print("EH = Energy (Hartrees), E = Energy (eV), NOrbGrp, OrbGrp, GrpDegen = degeneracies and corresponding orbital numbering by group in ePS, NormInt = single centre expansion convergence (should be ~1.0).")
    # #     print(molInfo['orbTable'].to_pandas()) # print() here doesn't give nice interactive result in notebook, although returning table does.
    #     orbPD = molInfo['orbTable'].to_pandas()
    #
    #     # Add symmery labels
    #     orbPD.insert(1, "Sym", molInfo['orbTable'].Sym)
    #
    #     # Tidy a bit...
    #     orbPD.drop(columns=["N","Type"], inplace=True)
    #
    #     # Check if notebook and output
    #     if isnotebook():
    #         display(orbPD) # Use IPython display(), works nicely for notebook output
    #         # Q: is display() always loaded automatically? I think so.
    #     else:
    #         print(orbPD)
    #
    #     # Check ExpOrb outputs...
    #     ind = (molInfo['orbTable'][:,8].values < 1-tolConv) + (molInfo['orbTable'][:,8].values > 1+tolConv)
    #     if ind.any():
    #         print(f"\n*** Warning: some orbital convergences outside single-center expansion convergence tolerance ({tolConv}):")
    #         print(molInfo['orbTable'][ind, [0, 8]].values)
    #
    # #     return orbPD  # Could also set return value here for simple orbPD printing.


    # def matEtoPD(self, keys = None, xDim = 'Eke', Erange = None, printTable = True, selDims = None, dType = None,
    #             thres = None, drop = True, fillna = False, squeeze = True, setPD = True):
    #     """
    #     Convert Xarray to PD for nice tabular display.
    #
    #     Basically code as per basicPlotters.lmPlot(), but looped over datasets.
    #
    #     """
    #
    #     # Default to all datasets
    #     if keys is None:
    #         keys = self.dataSets.keys()
    #
    #     pdConv = [] # List for outputs
    #
    #     for key in keys:
    #         # Init empty list for daPlotpd data
    #         if setPD:
    #             self.dataSets[key]['daPlotpd'] = []
    #
    #         for m, item in enumerate(self.dataSets[key]['matE']):
    #
    #             #*** Set, select & threshold
    #             da = self.dataSets[key]['matE'][m]
    #
    #             # Check xDim Eke/Ehv and clip if specified
    #             # More elegant way to swap on dims?
    #             if xDim == 'Ehv':
    #                 # Subset before plot to avoid errors on empty array selection!
    #                 da = da.swap_dims({'Eke':'Ehv'})
    #
    #             daSub = matEleSelector(da, thres=thres, inds = selDims, dims = xDim, drop = drop)
    #
    #             # NOTE: assumes xDim = Eke or Ehv, otherwise will fail.
    #             if Erange is not None:
    #                 # Check Etype, assuming only Eke or Ehv options.
    #                 Etype = 'Eke'
    #                 if 'Ehv' in daSub.dims:
    #                     Etype = 'Ehv'
    #
    #                 daSub = daSub.sel(**{Etype:slice(Erange[0], Erange[1])})  # With dict unpacking for var as keyword
    #
    #             #*** Data conversion if specified
    #             if dType is not None:
    #                 daSub = plotTypeSelector(daSub, pType = pType, axisUW = xDim)
    #
    #             #*** Convert to PD
    #             daPD, daSub = multiDimXrToPD(daSub, colDims = xDim, thres = thres, dropna = True, fillna = fillna, squeeze = squeeze)
    #
    #             pdConv.append(daPD)
    #
    #             if printTable:
    #                 print(f"\n*** {da.jobLabel}")
    #                 print(f"Matrix element table, threshold={thres}, data type={daSub.dtype}.")
    #                 # Check if notebook and output
    #                 if isnotebook():
    #                     display(daPD) # Use IPython display(), works nicely for notebook output
    #                     # Q: is display() always loaded automatically? I think so.
    #                 else:
    #                     print(daPD)
    #
    #             # Set Pandas table to dataset if specified.
    #             if setPD:
    #                 self.dataSets[key]['daPlotpd'].append(daPD)
    #
    #
    #     # Return value if not set to dataSets.
    #     if not setPD:
    #         return pdConv


# 14/10/20 Moved to base class
    # # **** Basic MFPAD plotting, see dev code in http://localhost:8888/lab/tree/dev/ePSproc/classDev/ePSproc_multijob_class_tests_N2O_011020_Stimpy.ipynb
    # def mfpadNumeric(self, selDims = {'Type':'L','it':1}, keys = None, res = 50):
    #     """
    #     MFPADs "direct" (numerical), without beta parameter computation.
    #
    #     Wrapper for :py:func:`epsproc.mfpad`, loops over all loaded datasets.
    #
    #     NOTE: for large datasets and/or large res, this can be memory-hungry.
    #
    #     """
    #
    #     # Default to all datasets
    #     if keys is None:
    #         keys = self.dataSets.keys()
    #
    #
    #     # Loop over datasets & store output
    #     for key in keys:
    #         TX = []
    #         for m, item in enumerate(self.dataSets[key]['matE']):
    #
    #             TXtemp, _ = mfpad(self.dataSets[key]['matE'][m], inds=selDims, res=res)  # Expand MFPADs
    #             TX.append(TXtemp)
    #
    #             # Propagate attrs
    #             TX[m].attrs = self.dataSets[key]['matE'][m].attrs
    #             TX[m]['dataType'] = 'mfpad'
    #
    #         self.dataSets[key]['TX'] = TX
    #
    #
    #
    # def mfpadPlot(self, selDims = {}, sumDims = {}, Erange = None, Etype = 'Eke', keys = None, pType = 'a', pStyle = 'polar', backend = 'mpl'):
    #
    #     # Default to all datasets
    #         if keys is None:
    #             keys = self.dataSets.keys()
    #         else:
    #             if isinstance(keys, int):   # Force list if int passed.
    #                 keys = [keys]
    #
    #         # Loop over datasets
    #         for key in keys:
    #             for m, item in enumerate(self.dataSets[key]['TX']):
    #                 plotFlag = True
    #
    #                 # More elegant way to swap on dims?
    #                 subset = self.dataSets[key]['TX'][m]
    #                 if Etype == 'Ehv':
    #                 # Subset before plot to avoid errors on empty array selection!
    #                     subset = subset.swap_dims({'Eke':'Ehv'})
    #
    #                 # Slice on Erange with/without step size - bit ugly.
    #                 if Erange is not None:
    #                     if len(Erange) == 3:
    #                         subset = subset.sel(**{Etype:slice(Erange[0], Erange[1], Erange[2])})   # With dict unpacking for var as keyword
    #                     else:
    #                         subset = subset.sel(**{Etype:slice(Erange[0], Erange[1])})
    #
    #                 # Handling for Euler or Labels dim
    #                 eDim = 'Euler'
    #                 if hasattr(subset, 'Labels'):
    # #                 if 'Labels' in selDims:
    #                     subset = subset.swap_dims({'Euler':'Labels'})
    #                     eDim = 'Labels'
    #
    #                 # if selDims is not None:  # Now set as empty dict to avoid issues later.
    #                 # TODO: smarter dim handling here - selDims may not be consistent over all dataSets!
    #                 if selDims:
    #                     subset = subset.sel(selDims)
    #
    #                 # Check dimensionality - sum over Sym & it if necessary
    #                 if subset.ndim > 3:
    #                     print(f"Found dims {subset.dims}, summing to reduce for plot. Pass selDims to avoid.")
    #                     if 'Sym' in subset.dims:
    #                         subset = subset.sum('Sym').squeeze()
    #                     if 'it' in subset.dims:
    #                         subset = subset.sum('it').squeeze()
    #
    #                     # Check result is OK
    #                     if (subset.ndim > 3) and not (eDim in subset.dims):
    #                         print(f"*** ERROR: dataset {self.dataSets[key]['TX'][m].jobLabel}  ndims = {subset.ndim}, dims = {subset.dims}. Skipping MFPAD plotting.")
    #                         plotFlag = False
    #
    #                     # else:
    #                     #     pass
    #
    #                 # Propagate attrs
    #                 subset.attrs = self.dataSets[key]['TX'][m].attrs
    #
    #                 # Plot
    #                 if plotFlag:
    #                     if eDim in subset.dims:
    #                         # TODO: make this more robust. Currently assumes Euler dim is 1st dim.
    #                         if subset[eDim].size > 1:
    #                             if pStyle is 'polar':
    #                                 for EulerInd in range(0,subset[eDim].size):
    #     #                                 print(f"*** Pol Geom: {subset[eDim][EulerInd].item()}")  # Need to pass all this through to plotters otherwise appears before plots!
    #
    #                                     _ = sphSumPlotX(subset[EulerInd], pType = pType, backend = backend, facetDim = Etype)
    #
    #                             elif pStyle is 'grid':
    #                                 print(f"Grid plot: {subset.attrs['jobLabel']}")
    #
    #                                 # Set data
    #                                 subset = plotTypeSelector(subset, pType = pType, axisUW = Etype)
    #
    #                                 if self.verbose['main']>2:
    #                                     print(subset)
    #
    #                                 # If Phi is sliced, assumed 2D (Theta,E) plots
    #                                 # NOTE - force real datatype here, otherwise get errors from complex residuals.
    #                                 # TODO: why is this (after abs pType selection?) Should error check & clean up - this could cause issues sometimes...
    #                                 if 'Phi' in selDims.keys():
    #                                     # print(subset)
    #                                     # subset.pipe(np.abs).plot(x='Theta', y=Etype, col=eDim, robust=True)
    #                                     subset.plot(x='Theta', y=Etype, col=eDim, robust=True)
    #                                 else:
    #                                     # subset.pipe(np.abs).plot(y='Theta',x='Phi', row=eDim, col=Etype, robust=True)
    #                                     subset.plot(y='Theta',x='Phi', row=eDim, col=Etype, robust=True)
    #
    #                                 # plt.title(subset.attrs['jobLabel'])
    #
    #                     else:
    #                         _ = sphSumPlotX(subset, pType = pType, backend = backend, facetDim = Etype)
