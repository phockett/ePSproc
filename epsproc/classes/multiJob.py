"""
Core classes for ePSproc data to handle multiple job filesets (energy and/or orbitals etc.)

23/09/20    Added ePSmultiJob class, currently very rough.

14/09/20    v1  Started class development. See ePSproc_multijob_class_dev_140920_bemo.ipynb

"""

from pathlib import Path
from matplotlib import pyplot as plt  # For plot legends with Xarray plotter
import xarray as xr
import pprint

# Local functions
from epsproc import readMatEle, headerFileParse, molInfoParse
from epsproc.util.summary import getOrbInfo

# Class for multiple ePS datasets
class ePSmultiJob():
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



    """

    # TODO: set to varg in for jobs dict
    def __init__(self, fileBase = None, jobDirs = None, jobStructure = None,
                 prefix = None, ext = '.out', Edp = 2, verbose = 1):

        self.verbose = verbose

        self.Edp = Edp  # Set default dp for Ehv conversion. May want to set this elsewhere...?

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

            if self.verbose > 0:
                print(f"Found ePS output files in subdirs: {testSub}")


        # Check root dir
        testRoot = list(self.jobs['fileBase'].glob(f"*{self.jobs['ext']}"))

        if testRoot and (self.jobs['jobStructure'] is None or 'rootDir'):
            self.jobs['jobDirs'] = [self.jobs['fileBase']]
            self.jobs['jobStructure'] = 'rootDir'

            if self.verbose > 0:
                print(f"Found ePS output files in root dir {self.jobs['fileBase']}.")


    def scanFiles(self):
        """
        Scan ePS output files from dir list.

        Adapted from https://phockett.github.io/ePSdata/XeF2-preliminary/XeF2_multi-orb_comparisons_270320-dist.html

        Currently outputting complicated dataSets dictionary/list - need to sort this out!
        - Entry per dir scanned.
        - Sub-entries per file, but collapsed in some cases.
        - Sending to Xarray dataset with orb labels should be cleaner.

        TODO:

        - convert outputs to Xarray dataset. Did this before, but currently missing file (on AntonJr)! CHECK BACKUPS - NOPE.
        - Confirm HV scaling - may be better to redo this, rather than correct existing values?

        """

        # Original version
        dataSets = {}  # Set dict to hold all data for the moment.
        jobs = {}  # Set dict to hold summary info

        # For Xarray Dataset stacking...
        dsXS = xr.Dataset()  # Set blank dataset
        dsMatE = xr.Dataset()


        for n, dirScan in enumerate(self.jobs['jobDirs']):
#             dataPath = Path(workingDir[0], dirScan)
            dataPath = dirScan  # Assume abs paths (?)

            # For dir scan
            dataSetXS = readMatEle(fileBase = dataPath, recordType = 'CrossSection', verbose = self.verbose)  # Set for XS + betas only
            dataSetMatE = readMatEle(fileBase = dataPath, recordType = 'DumpIdy', verbose = self.verbose)

            # Log some details - currently not passed directly from readMatEle()
            fN = len(dataSetXS)
            fList = [item.attrs['file'] for item in dataSetXS]

            if self.jobs['jobStructure'] == 'subDirs' and len(dataSetXS) > 1:
                if self.verbose > 2:
                    print('Processing data as subdir jobsets, with Eke stacking.')

                # Stack multi-E Xarrays into single array.
                # Keep results as list for compatibility with rest of code (otherwise will slice Xarray)
                dataXS = [xr.combine_nested(dataSetXS, concat_dim = ['Eke']).sortby('Eke')]
                dataMatE = [xr.combine_nested(dataSetMatE, concat_dim = ['Eke']).sortby('Eke')]

                # Set job info from first datafile of each set (orbital)
                dataFile = Path(dataXS[0].attrs['fileBase'], dataXS[0].attrs['file'])
                dataXS[0].attrs['jobInfo'] = headerFileParse(dataFile, verbose = self.verbose)
                dataXS[0].attrs['molInfo'] = molInfoParse(dataFile, verbose = self.verbose)

                # Set orb info
                dataXS[0].attrs['orbX'], dataXS[0].attrs['orbInfo'] = getOrbInfo(dataXS[0].attrs['jobInfo'], dataXS[0].attrs['molInfo'])

                # Set absolute photon energy
                dataXS[0]['Ehv'] = (dataXS[0]['Ehv'] - (float(dataXS[0].jobInfo['IPot']) + dataXS[0].orbX['E'].data[0])).round(self.Edp)
                dataMatE[0]['Ehv'] = (dataXS[0]['Ehv'] - (float(dataXS[0].jobInfo['IPot']) + dataXS[0].orbX['E'].data[0])).round(self.Edp)

                # Set as xr.ds(), label by orbital
                # TODO: fix this, in some cases gives errors with multiple values - probably an Ehv issue?
                # dsXS[dataXS[0].orbX['orb'].data[0]] = dataXS[0]
                # dsMatE[dataXS[0].orbX['orb'].data[0]] = dataMatE[0]

                # Job notes for plot labels etc.
                # This is basically as per summary.jobInfo() routine
                # TODO: Should replace/update/consolidate with regex and more to util func.
                if self.verbose > 1:
                    print(f"Batch: {dataXS[0].jobInfo['comments'][0].strip('#').strip()}")
                    print(f"Orbital: {dataXS[0].jobInfo['comments'][1].split(',', maxsplit=1)[1].strip()}")
                    print(f"OrbE: {dataXS[0].orbX['E'].data[0]}")

                jobNotes = []
                jobNotes.append({ 'batch': dataXS[0].jobInfo['comments'][0].strip('#').strip(),
                                'event': dataXS[0].jobInfo['comments'][1].split(',', maxsplit=1)[1].strip(),
                                'orbLabel': dataXS[0].jobInfo['comments'][1].split('(', maxsplit=1)[1].split(')')[0],
                                'orbE': dataXS[0].orbX['E'].data[0]
                                })

            else:
                if self.verbose > 2:
                    print('Processing data as single dir jobset, no Eke stacking.')

                dataXS = dataSetXS
                dataMatE = dataSetMatE
                jobNotes = []

                # Set job info for each file
                for m, item in enumerate(dataXS):

                    dataFile = Path(item.attrs['fileBase'], dataXS[m].attrs['file'])
                    dataXS[m].attrs['jobInfo'] = headerFileParse(dataFile, verbose = self.verbose)
                    dataXS[m].attrs['molInfo'] = molInfoParse(dataFile, verbose = self.verbose)

                    # Set orb info
                    dataXS[m].attrs['orbX'], dataXS[m].attrs['orbInfo'] = getOrbInfo(dataXS[m].attrs['jobInfo'], dataXS[m].attrs['molInfo'])

                    # Set absolute photon energy
                    dataXS[m]['Ehv'] = (dataXS[m]['Ehv'] - (float(dataXS[m].jobInfo['IPot']) + dataXS[m].orbX['E'].data[0])).round(self.Edp)
                    dataMatE[m]['Ehv'] = (dataXS[m]['Ehv'] - (float(dataXS[m].jobInfo['IPot']) + dataXS[m].orbX['E'].data[0])).round(self.Edp)

                    # Set as xr.ds(), label by orbital
                    dsXS[f"orb{dataXS[m].orbX['orb'].data[0]}"] = dataXS[m]
                    dsMatE[f"orb{dataXS[m].orbX['orb'].data[0]}"] = dataMatE[m]


                    # Job notes for plot labels etc.
                    # This is basically as per summary.jobInfo() routine
                    jobNotes.append({ 'batch': dataXS[m].jobInfo['comments'][0].strip('#').strip(),
                                    'event': dataXS[m].jobInfo['comments'][1].split(',', maxsplit=1)[1].strip(),
                                    'orbLabel': dataXS[m].jobInfo['comments'][1].split('(', maxsplit=1)[1].split(')')[0],
                                    'orbE': dataXS[m].orbX['E'].data[0]
                                    })

#                     print(m)
#                     print(dataXS[m].orbX['orb'].data[0])

#             # Set job info from first datafile of each set (orbital)
#             dataFile = Path(dataXS[0].attrs['fileBase'], dataXS[0].attrs['file'])
#             dataXS[0].attrs['jobInfo'] = ep.headerFileParse(dataFile)
#             dataXS[0].attrs['molInfo'] = ep.molInfoParse(dataFile)
#             dataXS.attrs['jobInfo'] = ep.headerFileParse(dataFile)
#             dataXS.attrs['molInfo'] = ep.molInfoParse(dataFile)

            dataSets[n] = {}
            dataSets[n]['dir'] = dirScan
            dataSets[n]['XS'] = dataXS
            dataSets[n]['matE'] = dataMatE
            dataSets[n]['jobNotes'] = jobNotes

            jobs[n] = {'dir': dirScan, 'fN': fN, 'fList': fList}

        self.dataSets = dataSets
        self.jobs['files'] = jobs

        # xr.ds types
        self.dsXS = dsXS
        self.dsMatE = dsMatE

        if self.verbose:
            self.jobsSummary()


    def jobsSummary(self):
        """
        Print some general info.

        TODO: add some info!
        """

        # Set pprint object
        pp = pprint.PrettyPrinter(indent=4)


        fN = [len(item['fList']) for item in self.jobs['files'].values()]
        print(f"Found {len(self.jobs['files'])} directories, with {sum(fN)} files.")
        print(f"Job structure: {self.jobs['jobStructure']}")

        if self.jobs['jobStructure'] == 'subDirs':
            print("(Results E stacked to one dataset per dir.)")

        # Print job details
#         [ep.jobSummary(testClass.dataSets[0]['XS'][1].jobInfo);
        for key in self.dataSets:
            print(f"\n*** Job dir {key} details")
            print(f"Directory: {self.jobs['files'][key]['dir']}")
            print(f"{self.jobs['files'][key]['fN']} files")

            for m, item in enumerate(self.dataSets[key]['XS']):
                # [print(line.strip('#')) for line in self.dataSets[key]['XS'][m].jobInfo['comments'][0:4]]
                # print(self.dataSets[key]['jobNotes'][m])
                pp.pprint(self.dataSets[key]['jobNotes'][m])




    # Define photon energy scale - now set at read-in.
#     def setHV(self):

#         # Set Ehv scale - currently incorrectly set by 1st vertical IP
#         for n, key in enumerate(self.dataSets):
#         #     dataSets[key][0]['Ehv'] -= (12.35 + orbIPs[n])
#             dataSets[key]['XS'][0]['Ehv'] = dataSets[key]['XS'][0]['Ehv'] - (12.35 + orbIPs[n])
#             dataSets[key]['matE'][0]['Ehv'] = dataSets[key]['matE'][0]['Ehv'] - (12.35 + orbIPs[n])

#             jobInfo['IPot']


    def plotGetCro(self, pType = 'SIGMA', Erange = None, Etype = 'Eke'):
        """
        Basic GetCro (cross-section) data plotting for multijob class. Run self.plot.line(x=Etype, col='Type') for each dataset.
        (See :py:func:`epsproc.classes.ePSmultiJob.plotGetCroComp()` for comparitive plots over datasets.)

        Note this is for LF averaged parameters, for more details see the `ePS starter notes <https://epsproc.readthedocs.io/en/latest/ePS_ePSproc_tutorial/ePS_tutorial_080520.html#Results>`_ for more details.

        Parameters
        ----------
        pType : str, optional, default = 'SIGMA'
            Set data for plotting, either 'SIGMA' (cross-section) or 'BETA' (B2 parameter).

        Erange : list of int or float, optional, default = None
            Set plot range [Emin, Emax]. Defaults to full data range if not set

        Etype : str, optional, default = 'Eke'
            Set plot dimension, either 'Eke' (electron kinetic energy) or 'Ehv' (photon energy).

        """

#         if self.jobs['jobStructure'] == 'subDirs':
        for key in self.dataSets:
#                 testClass.dataSets[key]['XS'][0].sel(XC='SIGMA', Eke=slice(Erange[0], Erange[1])).plot.line(x='Eke', col='Type')   # This works
            for m, item in enumerate(self.dataSets[key]['XS']):

                # Set default to full range, same for all cases
                if Erange is None:
                    Erange = [self.dataSets[key]['XS'][m][Etype].min().data, self.dataSets[key]['XS'][m][Etype].max().data]


                # More elegant way to swap on dims?
                if Etype == 'Ehv':
                # Subset before plot to avoid errors on empty array selection!
                    subset = self.dataSets[key]['XS'][m].swap_dims({'Eke':'Ehv'}).sel(XC=pType, **{Etype:slice(Erange[0], Erange[1])})   # With dict unpacking for var as keyword

                    # if subset.any():
                    #     subset.plot.line(x=Etype, col='Type')

                else:
                    subset = self.dataSets[key]['XS'][m].sel(XC=pType, **{Etype:slice(Erange[0], Erange[1])})   # With dict unpacking for var as keyword

                if subset.any():
                    subset.plot.line(x=Etype, col='Type')

                if pType == 'SIGMA':
                    plt.ylabel('XS/Mb')
                else:
                    plt.ylabel(r"$\beta_{LM}$")
#         else:

    # Basically as per plotGetCro, but subselect and put on single plot.
    # Should do all of this with Holoviews...!
    def plotGetCroComp(self, pType='SIGMA', pGauge='L', pSym=('All','All'), Erange = None, Etype = 'Eke'):
        """
        Basic GetCro (cross-section) data plotting for multijob class, comparitive plots.
        Run self.plot.line(x=Etype) for each dataset after subselection on Gauge and Symmetry, and use single axis.
        (See :py:func:`epsproc.classes.ePSmultiJob.plotGetCro()` for plots per dataset.)

        Note this is for LF averaged parameters, for more details see the `ePS starter notes <https://epsproc.readthedocs.io/en/latest/ePS_ePSproc_tutorial/ePS_tutorial_080520.html#Results>`_ for more details.

        Parameters
        ----------
        pType : str, optional, default = 'SIGMA'
            Set data for plotting, either 'SIGMA' (cross-section) or 'BETA' (B2 parameter).

        pGauge : str, optional, default = 'L'
            Set gauge, either 'L' (Length), 'V' (Velocity) or 'M' (Mixed)

        pSym : tuple of strs, optional, default = ('All','All')
            Select symmetry, (Cont, Targ).
            Default value will plot all allowed symmetries.

        Erange : list of int or float, optional, default = None
            Set plot range [Emin, Emax]. Defaults to full data range if not set

        Etype : str, optional, default = 'Eke'
            Set plot dimension, either 'Eke' (electron kinetic energy) or 'Ehv' (photon energy).

        """



        # Comparison plots over orbs

        # from matplotlib import pyplot as plt  # For legend
        lText = []


        for key in self.dataSets:
#                 testClass.dataSets[key]['XS'][0].sel(XC='SIGMA', Eke=slice(Erange[0], Erange[1])).plot.line(x='Eke', col='Type')   # This works
            for m, item in enumerate(self.dataSets[key]['XS']):

                # Set default to full range, same for all cases
                if Erange is None:
                    Erange = [self.dataSets[key]['XS'][m][Etype].min().data, self.dataSets[key]['XS'][m][Etype].max().data]


                # More elegant way to swap on dims?
                if Etype == 'Ehv':
                    # Subset before plot to avoid errors on empty array selection!
                    subset = self.dataSets[key]['XS'][m].swap_dims({'Eke':'Ehv'}).sel(XC=pType, Type=pGauge, Sym=pSym, **{Etype:slice(Erange[0], Erange[1])})  # With dict unpacking for var as keyword

                    # if subset.any():
                    #     pltObj = subset.plot.line(x=Etype)

                else:
                    subset = self.dataSets[key]['XS'][m].sel(XC=pType, Type=pGauge, Sym=pSym, **{Etype:slice(Erange[0], Erange[1])})   # With dict unpacking for var as keyword

                if subset.any():
                    pltObj = subset.plot.line(x=Etype)

                # Label with orb_sym
#                 lText.append(self.dataSets[key]['XS'][m].attrs['fileBase'].rsplit('/',maxsplit=1)[0])

                # lText.append(f"Orb {self.dataSets[key]['XS'][m].attrs['orbInfo']['orbN']} ({self.dataSets[key]['XS'][m].attrs['orbInfo']['orbSym'][0]})")

                lText.append(self.dataSets[key]['jobNotes'][m]['orbLabel'])

        plt.legend(lText)

        if pType == 'SIGMA':
            plt.ylabel('XS/Mb')
        else:
            plt.ylabel(r"$\beta_{LM}$")


    def lmPlot(self, Erange = None, Etype = 'Eke', refData = None, **kwargs):
        """
        Wrapper for :py:func:`epsproc.lmPlot` for multijob class. Run lmPlot() for each dataset.

        Parameters
        ----------
        Erange : list of int or float, optional, default = None
            Set plot range [Emin, Emax]. Defaults to full data range if not set

        Etype : str, optional, default = 'Eke'
            Set plot dimension, either 'Eke' (electron kinetic energy) or 'Ehv' (photon energy).

        refData : tuple (key,m), optional, default = None
            If set, calculate difference plots against reference dataset.
            TODO: implement difference plots.
            TODO: implement testing logic, may fail without E-axis forcing, and sym summation?

        kwargs : dict, optional, default = {}
            Plotting options to pass to :py:func:`epsproc.lmPlot`.
            These will also be set in self.lmPlotOpts for further use.
            Note that any existing options in self.lmPlotOpts will also be used, or overwritten if matching keys are found.


        Notes
        -----
        Basic scheme from ePSmultijob.plotGetCro, which loops and switches on Eke/Ehv. Should tidy up at some point.

        """

        # Set lmPlotOpts
        # Check passed args vs. self.lmPlotOpts and overwrite
        if kwargs:
            for key, value in kwargs.items():
                self.lmPlotOpts[key] = value


        # Loop over datasets
        for key in self.dataSets:
    #                 testClass.dataSets[key]['XS'][0].sel(XC='SIGMA', Eke=slice(Erange[0], Erange[1])).plot.line(x='Eke', col='Type')   # This works

            # Init empty list for daPlotpd data
            self.dataSets[key]['daPlotpd'] = []

            for m, item in enumerate(self.dataSets[key]['matE']):

                # Set default to full range, same for all cases
                if Erange is None:
                    Erange = [self.dataSets[key]['matE'][m][Etype].min().data, self.dataSets[key]['matE'][m][Etype].max().data]


                # More elegant way to swap on dims?
                if Etype == 'Ehv':
                # Subset before plot to avoid errors on empty array selection!
                    subset = self.dataSets[key]['matE'][m].swap_dims({'Eke':'Ehv'}).sel(**{Etype:slice(Erange[0], Erange[1])})   # With dict unpacking for var as keyword

                else:
                    subset = self.dataSets[key]['matE'][m].sel(**{Etype:slice(Erange[0], Erange[1])})   # With dict unpacking for var as keyword

                # Run lmPlot
                if subset.any():
                    daPlot, daPlotpd, legendList, gFig = ep.lmPlot(subset, xDim = Etype, **self.lmPlotOpts)
                else:
                    daPlotpd = None

                self.dataSets[key]['daPlotpd'].append(daPlotpd)  # Set to include None cases to keep indexing. Should set as dict instead?
