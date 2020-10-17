"""
Core classes for ePSproc data.

13/10/20    v1  Started class development, reverse engineering a little from multiJob.py case.


"""
import pprint
# import os
from pathlib import Path
from matplotlib import pyplot as plt  # For plot legends with Xarray plotter
import numpy as np # Needed only for np.nan at the moment

# Local functions
from epsproc import readMatEle, headerFileParse, molInfoParse, lmPlot, matEleSelector, plotTypeSelector, multiDimXrToPD, mfpad, sphSumPlotX
# from epsproc.classes.base import ePSbase
from epsproc.util.summary import getOrbInfo, molPlot
from epsproc.util.env import isnotebook

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

        # Set for main functions and subfunctions
        self.verbose = {'main':verbose, 'sub':verbose-1}

        if self.verbose['sub'] < 0:
            self.verbose['sub'] = 0


        self.Edp = Edp  # Set default dp for Ehv conversion. May want to set this elsewhere, and match to data resolution?
        self.lmPlotOpts = {}  # Set empty dict for plot options.

        # If not set, try working dir
        if fileBase is None:
#             fileBase = os.getcwd()
            fileBase = Path.cwd()

        # Set file properties
        self.job = {'fileBase':Path(fileBase),
                     'ext':ext,
                     'prefix':prefix,
                     }

        # Set master dictionary to hold sorted datasets
        # Set here to allow persistence over multiple scanFiles() runs from child classes
        self.data = {}
        self.jobNotes = []
#         self.jobs = []  # Possibly redundant, but usful for multi dir case


    def scanFiles(self, dataPath = None, reset = False, keyType = 'orb'):
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

        """

        # Allow override here for subclasses/independent use
        if dataPath is None:
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

        # Read fileset
        # 13/10/20 with updated sorting code, this should return
        # - a one-element list for a dir with Eke split files.
        # - a multi-element list for a dir with multiple jobs.
        # - Note cross-over with multiJob class in latter case.
        dataSetXS = readMatEle(fileBase = dataPath, recordType = 'CrossSection', verbose = self.verbose['sub'])  # Set for XS + betas only
        dataSetMatE = readMatEle(fileBase = dataPath, recordType = 'DumpIdy', verbose = self.verbose['sub'])

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

            # Additional labels, use these in plotting routines later
            dataSetXS[m].attrs['jobLabel'] = item.jobInfo['comments'][1].split('(', maxsplit=1)[1].split(')')[0]
            dataSetMatE[m].attrs['jobLabel'] = item.jobInfo['comments'][1].split('(', maxsplit=1)[1].split(')')[0]

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
            else:
                key = m

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

            dsSet[key]['job'] = {'dir': dataPath, 'fN': fN, 'files': fList}
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
            print(f"Dir {self.data[key]['job']['dir']}, {self.data[key]['job']['fN']} files.")
            pp.pprint(self.data[key]['jobNotes'])
#             print(f"Directory: {self.jobs['files'][key]['dir']}")
#             print(f"{self.job['files']['fN']} files")

            if self.verbose['main'] > 1:
                print("File list: ", *self.data[key]['job']['files'], sep='\n   ')

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
        orbPD.insert(1, "Sym", molInfo['orbTable'].Sym)

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

    # **** BLM calculation wrappers

    def AFBLM(self, keys, **kwargs):
        """
        Basic wrapper for :py:func:`epsproc.geomFunc.afblmXprod`.

        Currently set to calculate for all data, with afblmXprod defaults, unless additional kwargs passed.

        TODO:
        - Add subselection here?

        """

        # Default to all datasets
        if keys is None:
            keys = list(self.data.keys())

        for key in keys:
            if self.verbose['main']:
                print(f"\nCalculating AF-BLMs for job key: {key}")

            _, _, _, BetaNormX = ep.geomFunc.afblmXprod(self.data[key]['matE'], **kwargs)

            self.data[key]['AFBLM'] = BetaNormX


    # ************** Plotters

    def plotGetCro(self, pType = 'SIGMA', Erange = None, Etype = 'Eke', keys = None, backend = 'mpl'):
        """
        Basic GetCro (cross-section) data plotting for multijob class. Run self.plot.line(x=Etype, col='Type') for each dataset.
        (See :py:func:`epsproc.classes.ePSmultiJob.plotGetCroComp()` for comparitive plots over datasets.)

        Note this is for LF averaged parameters, for more details see the `ePS starter notes <https://epsproc.readthedocs.io/en/latest/ePS_ePSproc_tutorial/ePS_tutorial_080520.html#Results>`_ for more details.

        Parameters
        ----------
        pType : str, optional, default = 'SIGMA'
            Set data for plotting, either 'SIGMA' (cross-section) or 'BETA' (B2 parameter).
            If backend = 'hv' this parameter is not used.

        Erange : list of int or float, optional, default = None
            Set plot range [Emin, Emax]. Defaults to full data range if not set

        Etype : str, optional, default = 'Eke'
            Set plot dimension, either 'Eke' (electron kinetic energy) or 'Ehv' (photon energy).

        keys : list, optional, default = None
            Keys for datasets to plot.
            If None, all datasets will be plotted.

        backend : str, optional, default = 'mpl'
            Set plotter to use.

            - 'mpl' : Use Matplotlib/native Xarray plotter
            - 'hv' : use Holoviews via :py:func:`epsproc.plotters.hvPlotters.XCplot()`

        """
        # Default to all datasets
        if keys is None:
            keys = self.data.keys()

#         if self.jobs['jobStructure'] == 'subDirs':
        for key in keys:
#                 testClass.dataSets[key]['XS'][0].sel(XC='SIGMA', Eke=slice(Erange[0], Erange[1])).plot.line(x='Eke', col='Type')   # This works
            # for m, item in enumerate(self.data[key]['XS']):

            # Set default to full range, same for all cases
            if Erange is None:
                Erange = [self.data[key]['XS'][Etype].min().data, self.data[key]['XS'][Etype].max().data]

            subset = self.data[key]['XS']
            jobLabel = self.data[key]['XS'].jobLabel

            # More elegant way to swap on dims?
            if Etype == 'Ehv':
            # Subset before plot to avoid errors on empty array selection!
                subset = subset.swap_dims({'Eke':'Ehv'}) #.sel(**{Etype:slice(Erange[0], Erange[1])})   # With dict unpacking for var as keyword

                # if subset.any():
                #     subset.plot.line(x=Etype, col='Type')

            subset = subset.sel(**{Etype:slice(Erange[0], Erange[1])})   # With dict unpacking for var as keyword

            if subset.any():
                # Propagate attrs for plot labelling
                subset.attrs = self.data[key]['XS'].attrs

                if backend == 'mpl':
                    subset.sel(XC=pType).plot.line(x=Etype, col='Type')

                    if pType == 'SIGMA':
                        plt.ylabel('XS/Mb')
                    else:
                        plt.ylabel(r"$\beta_{LM}$")

                if backend == 'hv':
                    layout, *_ = hvPlotters.XCplot(subset, kdims = Etype)

                    print(f"\n*** {jobLabel}")
                    # Check if notebook and output
                    if isnotebook():
                        display(layout) # Use IPython display(), works nicely for notebook output
                        # Q: is display() always loaded automatically? I think so.
                    else:
                        layout  # Not yet tested - not sure what the options are here for hv/Bokeh.
                                # May just want to return this object (hv layout)?

                    # Setting title here doesn't work, now set in XCplot()
                    # display(layout.opts(title=jobLabel))

#         else:


    # Basically as per plotGetCro, but subselect and put on single plot.
    # Should do all of this with Holoviews...!
    def plotGetCroComp(self, pType='SIGMA', pGauge='L', pSym=('All','All'), Erange = None, Etype = 'Eke',  keys = None):
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

        keys : list, optional, default = None
            Keys for datasets to plot.
            If None, all datasets will be plotted.

        """



        # Comparison plots over orbs

        # from matplotlib import pyplot as plt  # For legend
        lText = []

        # Default to all datasets
        if keys is None:
            keys = self.data.keys()

        for key in keys:
#                 testClass.dataSets[key]['XS'][0].sel(XC='SIGMA', Eke=slice(Erange[0], Erange[1])).plot.line(x='Eke', col='Type')   # This works
            # for m, item in enumerate(self.data[key]['XS']):

            # Set default to full range, same for all cases
            if Erange is None:
                Erange = [self.data[key]['XS'][Etype].min().data, self.data[key]['XS'][Etype].max().data]


            # More elegant way to swap on dims?
            if Etype == 'Ehv':
                # Subset before plot to avoid errors on empty array selection!
                subset = self.data[key]['XS'].swap_dims({'Eke':'Ehv'}).sel(XC=pType, Type=pGauge, Sym=pSym, **{Etype:slice(Erange[0], Erange[1])})  # With dict unpacking for var as keyword

                # if subset.any():
                #     pltObj = subset.plot.line(x=Etype)

            else:
                subset = self.data[key]['XS'].sel(XC=pType, Type=pGauge, Sym=pSym, **{Etype:slice(Erange[0], Erange[1])})   # With dict unpacking for var as keyword

            if subset.any():
                pltObj = subset.plot.line(x=Etype)
                lText.append(self.data[key]['jobNotes']['orbLabel'])

            # Label with orb_sym
#                 lText.append(self.dataSets[key]['XS'][m].attrs['fileBase'].rsplit('/',maxsplit=1)[0])

            # lText.append(f"Orb {self.dataSets[key]['XS'][m].attrs['orbInfo']['orbN']} ({self.dataSets[key]['XS'][m].attrs['orbInfo']['orbSym'][0]})")

            # lText.append(self.dataSets[key]['jobNotes'][m]['orbLabel'])

        plt.legend(lText)

        if pType == 'SIGMA':
            plt.ylabel('XS/Mb')
        else:
            plt.ylabel(r"$\beta_{LM}$")


    def lmPlot(self, Erange = None, Etype = 'Eke', dataType = 'matE', xDim = None, keys = None, refDataKey = None, reindexTol = 0.5, reindexFill = np.nan, setPD = True, **kwargs):
        """
        Wrapper for :py:func:`epsproc.lmPlot` for multijob class. Run lmPlot() for each dataset.

        Parameters
        ----------
        Erange : list of int or float, optional, default = None
            Set plot range [Emin, Emax]. Defaults to full data range if not set

        Etype : str, optional, default = 'Eke'
            Set plot dimension, either 'Eke' (electron kinetic energy) or 'Ehv' (photon energy).

        dataType : str, optional, default = 'matE'
            Set data type to plot, corresponding to label in self.data
            - 'matE' raw matrix elements.
            - 'AFBLM' computed AF BLMs.

        xDim : str, optional, default = None
            Settings for x-axis, if None plot vs. Etype.
            See :py:func:`epsproc.lmPlot` for more details.

        keys : list, optional, default = None
            Keys for datasets to plot.
            If None, all datasets will be plotted.

        refDataKey : str or int, optional, default = None
            If set, calculate difference plots against reference dataset.
            This must be a key in self.data.
            TODO: implement difference plots.
            TODO: implement testing logic, may fail without E-axis forcing, and sym summation?

        reindexTol : float, optional, default = 0.1
            If computing difference data, the reference data is reindexed to ensure E grid matching.
            This specifies tolerance (in E units, usually eV) for reindexing.
            If this fails, difference plot may be null.

        reindexFill : int or float, optional, default = NaN
            Value to use for missing values upon reindexing.
            Default matches [Xarray.reindex default](http://xarray.pydata.org/en/stable/generated/xarray.DataArray.reindex.html), i.e. NaN, but this may give issues in some cases.

        setPD : bool, optional, default = True
            Set Pandas array in main dataset?

        kwargs : dict, optional, default = {}
            Plotting options to pass to :py:func:`epsproc.lmPlot`.
            These will also be set in self.lmPlotOpts for further use.
            Note that any existing options in self.lmPlotOpts will also be used, or overwritten if matching keys are found.


        Notes
        -----
        Basic scheme from ePSmultijob.plotGetCro, which loops and switches on Eke/Ehv. Should tidy up at some point.

        """

        # Default to all datasets
        if keys is None:
            keys = list(self.data.keys())

        # Set lmPlotOpts
        # Check passed args vs. self.lmPlotOpts and overwrite
        if kwargs:
            for key, value in kwargs.items():
                self.lmPlotOpts[key] = value

        # Set default to full range of 1st dataset, keep same for all cases
        # TODO: set per dataset?
        if Erange is None:
            Erange = [self.data[keys[0]][dataType][Etype].min().data, self.data[keys[0]][dataType][Etype].max().data]

        # Set ref dataset if required
        if refDataKey is not None:
            refData = self.data[refDataKey][dataType]

            if Etype == 'Ehv':
                refData = refData.swap_dims({'Eke':'Ehv'})

            refData = refData.sel(**{Etype:slice(Erange[0], Erange[1])})
            refData.attrs = self.data[refDataKey][dataType].attrs # Propagate atrrs.

        else:
            refData = None


        # Loop over datasets
        for key in keys:
    #                 testClass.dataSets[key]['XS'][0].sel(XC='SIGMA', Eke=slice(Erange[0], Erange[1])).plot.line(x='Eke', col='Type')   # This works

            # Init empty list for daPlotpd data
            if setPD:
                self.data[key]['daPlotpd'] = []

            # for m, item in enumerate(self.data[key]['matE']):

            # More elegant way to swap on dims?
            if Etype == 'Ehv':
            # Subset before plot to avoid errors on empty array selection!
                subset = self.data[key][dataType].swap_dims({'Eke':'Ehv'}).sel(**{Etype:slice(Erange[0], Erange[1])})   # With dict unpacking for var as keyword

            else:
                subset = self.data[key][dataType].sel(**{Etype:slice(Erange[0], Erange[1])})   # With dict unpacking for var as keyword

            # Difference data
            # NOTE: ref. data is reindexed here to ensure E-point subtraction (otherwise mismatch will give null results)
            # TODO: may want some selection logic here, this may produce empty results in many cases.
            if refData is not None:
                subset = subset - refData.reindex(**{Etype:subset[Etype]}, method='nearest', tolerance=reindexTol, fill_value=reindexFill)
                subset.attrs = self.data[key][dataType].attrs  # Propagate attribs
                subset.attrs['jobLabel'] = f"{subset.attrs['jobLabel']} diff with {refData.attrs['jobLabel']}"

                # Slightly crude check for empty result.
                # Note this is different from subset.any() below, which only checks for 0
                # See https://numpy.org/doc/stable/reference/generated/numpy.any.html
                if subset.max().isnull():
                    print("*** Warning, difference array is Null. Try a larger reindexTol value.")

            # Run lmPlot
            if subset.any():
                daPlot, daPlotpd, legendList, gFig = lmPlot(subset, xDim = xDim, **self.lmPlotOpts)
            else:
                daPlotpd = None

            # Set Pandas table to dataset if specified.
            if setPD:
                self.data[key]['daPlotpd'].append(daPlotpd)  # Set to include None cases to keep indexing. Should set as dict instead?



    # **** Basic MFPAD numerics & plotting, see dev code in http://localhost:8888/lab/tree/dev/ePSproc/classDev/ePSproc_multijob_class_tests_N2O_011020_Stimpy.ipynb
    def mfpadNumeric(self, selDims = {'Type':'L','it':1}, keys = None, res = 50):
        """
        MFPADs "direct" (numerical), without beta parameter computation.

        Wrapper for :py:func:`epsproc.mfpad`, loops over all loaded datasets.

        NOTE: for large datasets and/or large res, this can be memory-hungry.

        """

        # Default to all datasets
        if keys is None:
            keys = self.data.keys()


        # Loop over datasets & store output
        for key in keys:
            # TX = []
            # for m, item in enumerate(self.dataSets[key]['matE']):

            TX, _ = mfpad(self.data[key]['matE'], inds=selDims, res=res)  # Expand MFPADs
            # TX.append(TXtemp)

            # Propagate attrs
            TX.attrs = self.data[key]['matE'].attrs
            TX['dataType'] = 'mfpad'

            self.data[key]['TX'] = TX



    def mfpadPlot(self, selDims = {}, sumDims = {}, Erange = None, Etype = 'Eke', keys = None, pType = 'a2', pStyle = 'polar', backend = 'mpl'):

        # Default to all datasets
            if keys is None:
                keys = self.data.keys()
            else:
                if not isinstance(keys, list):   # Force list if single item passed
                    keys = [keys]

            # Loop over datasets
            for key in keys:
                # for m, item in enumerate(self.dataSets[key]['TX']):
                plotFlag = True

                # More elegant way to swap on dims?
                subset = self.data[key]['TX']
                if Etype == 'Ehv':
                # Subset before plot to avoid errors on empty array selection!
                    subset = subset.swap_dims({'Eke':'Ehv'})

                # Slice on Erange with/without step size - bit ugly.
                if Erange is not None:
                    if len(Erange) == 3:
                        subset = subset.sel(**{Etype:slice(Erange[0], Erange[1], Erange[2])})   # With dict unpacking for var as keyword
                    else:
                        subset = subset.sel(**{Etype:slice(Erange[0], Erange[1])})

                # Handling for Euler or Labels dim
                eDim = 'Euler'
                if hasattr(subset, 'Labels'):
#                 if 'Labels' in selDims:
                    subset = subset.swap_dims({'Euler':'Labels'})
                    eDim = 'Labels'

                # if selDims is not None:  # Now set as empty dict to avoid issues later.
                # TODO: smarter dim handling here - selDims may not be consistent over all dataSets!
                if selDims:
                    subset = subset.sel(selDims)

                # Check dimensionality - sum over Sym & it if necessary
                if subset.ndim > 3:
                    print(f"Found dims {subset.dims}, summing to reduce for plot. Pass selDims to avoid.")
                    if 'Sym' in subset.dims:
                        subset = subset.sum('Sym').squeeze()
                    if 'it' in subset.dims:
                        subset = subset.sum('it').squeeze()

                    # Check result is OK
                    if (subset.ndim > 3) and not (eDim in subset.dims):
                        print(f"*** ERROR: dataset {self.data[key]['TX'].jobLabel}  ndims = {subset.ndim}, dims = {subset.dims}. Skipping MFPAD plotting.")
                        plotFlag = False

                    # else:
                    #     pass

                # Propagate attrs
                subset.attrs = self.data[key]['TX'].attrs

                # Plot
                if plotFlag:
                    if eDim in subset.dims:
                        # TODO: make this more robust. Currently assumes Euler dim is 1st dim.
                        if subset[eDim].size > 1:
                            if pStyle is 'polar':
                                for EulerInd in range(0,subset[eDim].size):
    #                                 print(f"*** Pol Geom: {subset[eDim][EulerInd].item()}")  # Need to pass all this through to plotters otherwise appears before plots!

                                    _ = sphSumPlotX(subset[EulerInd], pType = pType, backend = backend, facetDim = Etype)

                            elif pStyle is 'grid':
                                print(f"Grid plot: {subset.attrs['jobLabel']}")

                                # Set data
                                subset = plotTypeSelector(subset, pType = pType, axisUW = Etype)

                                if self.verbose['main']>2:
                                    print(subset)

                                # If Phi is sliced, assumed 2D (Theta,E) plots
                                # NOTE - force real datatype here, otherwise get errors from complex residuals.
                                # TODO: why is this (after abs pType selection?) Should error check & clean up - this could cause issues sometimes...
                                if 'Phi' in selDims.keys():
                                    # print(subset)
                                    # subset.pipe(np.abs).plot(x='Theta', y=Etype, col=eDim, robust=True)
                                    subset.plot(x='Theta', y=Etype, col=eDim, robust=True)
                                else:
                                    # subset.pipe(np.abs).plot(y='Theta',x='Phi', row=eDim, col=Etype, robust=True)
                                    subset.plot(y='Theta',x='Phi', row=eDim, col=Etype, robust=True)

                                # plt.title(subset.attrs['jobLabel'])

                    else:
                        _ = sphSumPlotX(subset, pType = pType, backend = backend, facetDim = Etype)
