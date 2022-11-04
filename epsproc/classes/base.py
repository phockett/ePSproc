"""
Core classes for ePSproc data.

12/10/22    Added phaseCorrection() wrapper.

13/10/20    v1  Started class development, reverse engineering a little from multiJob.py case.


TODO:

- Centralise subselection function.
- Error checking on datatype per key for plotters etc.

"""
import pprint
# import os
from pathlib import Path
from matplotlib import pyplot as plt  # For plot legends with Xarray plotter
import numpy as np # Needed only for np.nan at the moment
import scipy.constants

# Local functions
# from epsproc import lmPlot, matEleSelector, plotTypeSelector, multiDimXrToPD, mfpad, sphSumPlotX, sphFromBLMPlot
# from epsproc import readMatEle, headerFileParse, molInfoParse, lmPlot, matEleSelector, plotTypeSelector, multiDimXrToPD, mfpad, sphSumPlotX, sphFromBLMPlot
# from epsproc.classes.base import ePSbase
# from epsproc.util.summary import getOrbInfo, molPlot
# from epsproc.util.env import isnotebook

from epsproc import mfpad, plotTypeSelector, writeXarray, readXarray
from epsproc.MFPAD import mfWignerDelay
from epsproc.geomFunc.afblmGeom import afblmXprod, AFwfExp
from epsproc.geomFunc.mfblmGeom import mfblmXprod
from epsproc.calc.phases import phaseCorrection

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
    - Stick with .data for all data, or just promote to top-level keys per dataset? This might be neater overall.
    - Change to decorators for basic function wrappers - should be cleaner, and enforce method/style.
    - Check file IO logic, some of this is already handled in lower level codes.


    """

    # Import methods - these are essentially wrappers for core functions
    from ._IO import scanFiles, jobsSummary, molSummary, matEtoPD
    from ._plotters import BLMplot, padPlot, lmPlot, plotGetCro, plotGetCroComp, ADMplot, _hvBLMplot
    from ._selectors import Esubset

    # TODO: set to varg in for jobs dict
    def __init__(self, fileBase = None, fileIn = None,
                 prefix = None, ext = '.out', Edp = 1, verbose = 1,
                 thres = 1e-2, thresDims = 'Eke', selDims = {'Type':'L'}):

        # Set for main functions and subfunctions
        self.verbose = {'main':verbose, 'sub':verbose-1}

        if self.verbose['sub'] < 0:
            self.verbose['sub'] = 0

        # Calculation options
        self.Edp = Edp  # Set default dp for Ehv conversion. May want to set this elsewhere, and match to data resolution?
        self.calcOpts = {'thres':thres, 'thresDims':thresDims, 'selDims':selDims}

        # Plotting options
        self.lmPlotOpts = {'SFflag':False}  # Set empty dict for plot options.

        # If not set, try working dir
        if fileBase is None:
#             fileBase = os.getcwd()
            fileBase = Path.cwd()
        else:
            fileBase = Path(fileBase)  # Force Path object

        # if fileIn is not None:
        #     fileIn = Path(fileIn)

        # Set file properties
        self.job = {'fileBase':fileBase,
                    'fileIn':fileIn,
                     'ext':ext,
                     'prefix':prefix,
                     }

        # Set master dictionary to hold sorted datasets
        # Set here to allow persistence over multiple scanFiles() runs from child classes
        self.data = {}
        self.jobNotes = []
#         self.jobs = []  # Possibly redundant, but usful for multi dir case

        # 04/11/22 - add plot bag.
        # May want more flexibility here?
        self.plots = {}


    # **** Small utility fns.
    def _keysCheck(self, keys):
        """
        Set keys:

        - Default to all if None.
        - Force list if singleton key passed.

        """
        # Default to all datasets
        if keys is None:
            keys = list(self.data.keys())  # Get keys and set as list
        else:
            if not isinstance(keys, list):   # Force list if single item passed
                keys = [keys]

        return keys

    # **** BLM calculation wrappers
    def MFBLM(self, keys = None, **kwargs):
        """
        Basic wrapper for :py:func:`epsproc.geomFunc.mfblmXprod`.

        Currently set to calculate for all data, with mfblmXprod defaults, unless additional kwargs passed.

        TODO:
        - Add subselection here?

        """

        # # Default to all datasets
        # if keys is None:
        #     keys = list(self.data.keys())
        keys = self._keysCheck(keys)

        for key in keys:
            if self.verbose['main']:
                print(f"\nCalculating MF-BLMs for job key: {key}")

            BetaNormX = mfblmXprod(self.data[key]['matE'], **kwargs)

            self.data[key]['MFBLM'] = BetaNormX




    def AFBLM(self, keys = None, **kwargs):
        """
        Basic wrapper for :py:func:`epsproc.geomFunc.afblmXprod`.

        Currently set to calculate for all data, with afblmXprod defaults, unless additional kwargs passed.

        TODO:
        - Add subselection here?

        """

        # # Default to all datasets
        # if keys is None:
        #     keys = list(self.data.keys())
        keys = self._keysCheck(keys)

        for key in keys:
            if self.verbose['main']:
                print(f"\nCalculating AF-BLMs for job key: {key}")

            BetaNormX = afblmXprod(self.data[key]['matE'], **kwargs)

            self.data[key]['AFBLM'] = BetaNormX


    # **** Basic AFPAD numerics, see dev code in http://localhost:8888/lab/tree/SLAC/angular_streaking/AF_wavefns_method_dev_050121.ipynb
    # Function here basically as per mfpadNumeric
    def afpadNumeric(self, keys = None, **kwargs):
        """
        AFPADs "direct" (numerical), without beta parameter computation.

        Wrapper for :py:func:`epsproc.afblmGeom.AFwfExp()`, loops over all loaded datasets.

        NOTE: this is preliminary and unverified.

        Parameters
        ----------
        keys : str, int or list, optional, default = None
            If set, use only these datasets (keys).
            Otherwise run for all datasets in self.data.

        **kwargs : optional
            Args passed to :py:func:`epsproc.mfpad`.

        NOTE: for large datasets and/or large res, this can be memory-hungry.

        Returns
        -------

        None, but sets self.data[key]['TXaf'] and self.data[key]['DeltaKQS'] for each key.

        """

        # # Default to all datasets
        # if keys is None:
        #     keys = self.data.keys()
        keys = self._keysCheck(keys)


        # Loop over datasets & store output
        for key in keys:
            # TX = []
            # for m, item in enumerate(self.dataSets[key]['matE']):

            TX, AFterm, DeltaKQS = AFwfExp(self.data[key]['matE'], **kwargs)  # Expand MFPADs
            # TX.append(TXtemp)

            # Propagate attrs
            TX.attrs = self.data[key]['matE'].attrs
            TX.attrs['dataType'] = 'afpad'

            DeltaKQS.attrs = TX.attrs
            DeltaKQS.attrs['dataType'] = 'DeltaKQS'

            self.data[key]['TXaf'] = TX
            self.data[key]['DeltaKQS'] = DeltaKQS

            if self.verbose['main']:
                print(f'Set AF expansion parameters TXaf and DeltaKQS for {key}')




    # **** Basic MFPAD numerics & plotting, see dev code in http://localhost:8888/lab/tree/dev/ePSproc/classDev/ePSproc_multijob_class_tests_N2O_011020_Stimpy.ipynb
    def mfpadNumeric(self, keys = None, **kwargs):
        """
        MFPADs "direct" (numerical), without beta parameter computation.

        Wrapper for :py:func:`epsproc.mfpad`, loops over all loaded datasets.

        Parameters
        ----------
        keys : str, int or list, optional, default = None
            If set, use only these datasets (keys).
            Otherwise run for all datasets in self.data.

        **kwargs : optional
            Args passed to :py:func:`epsproc.mfpad`.

        NOTE: for large datasets and/or large res, this can be memory-hungry.

        """

        # # Default to all datasets
        # if keys is None:
        #     keys = self.data.keys()
        keys = self._keysCheck(keys)


        # Loop over datasets & store output
        for key in keys:
            # TX = []
            # for m, item in enumerate(self.dataSets[key]['matE']):

            TX, _ = mfpad(self.data[key]['matE'], **kwargs)  # Expand MFPADs
            # TX.append(TXtemp)

            # Propagate attrs
            TX.attrs = self.data[key]['matE'].attrs
            TX.attrs['dataType'] = 'mfpad'

            self.data[key]['TX'] = TX



    def wignerDelay(self, keys = None, pType = 'phaseUW', selDims={'Type':'L','it':1}, **kwargs):
        """
        Wigner delay computation as phase derivative of TX grid.

        Multi data-set wrapper for numerics; uses :py:func:`epsproc.MFPAD.mfpad()` and :py:func:`epsproc.MFPAD.mfWignerDelay()`.

        27/10/20 initial version added. Looks OK for N2 & N2O test cases, but not carefully tested as yet.
        http://localhost:8888/lab/tree/dev/ePSproc/classDev/ePSproc_class_demo_161020_Wigner_271020.ipynb


        Parameters
        ----------
        keys : str, int or list, optional, default = None
            If set, use only these datasets (keys).
            Otherwise run for all datasets in self.data.

        pType : str, optional, default = 'phaseUW'
            Used to set data conversion options, as implemented in :py:func:`epsproc.plotTypeSelector()`
            - 'phase' use np.angle()
            - 'phaseUW' unwapped with np.unwrap(np.angle())

        **kwargs : optional
            Args passed to :py:func:`epsproc.mfpad`.


        """

        keys = self._keysCheck(keys)

        # Compute energy derivatives
        for key in keys:

            # Set full wavefn expansion if not present.
            if 'TX' not in self.data[key].keys():
                # self.mfpadNumeric(selDims = selDims, keys = key, **kwargs)
                self.mfpadNumeric(keys = key, inds = selDims, **kwargs)   # 07/07/22 - need to rationalise selDims here! Now matched basic MFPAD.mfpad() routine

            # self.data[key]['wigner'] = self.data[key]['TX'].copy() # Init array with copy (may not be necessary)
            #
            # self.data[key]['wigner']['EJ'] = self.data[key]['wigner']['Eke']*scipy.constants.e  # Convert to J units
            #
            # unitConv = scipy.constants.hbar/1e-18  # Convert to attoseconds
            # # TODO: move this to util functions.
            #
            # # Use existing functionality to set phases - easier if unwrap required (uses lambdas)
            # self.data[key]['wigner'] = plotTypeSelector(self.data[key]['wigner'].conj(), pType = 'phaseUW').differentiate('EJ')* unitConv # Same results as above, bit smoother for phaseUW case
            #
            # # Propagate attrs
            # self.data[key]['wigner'].attrs = self.data[key]['TX'].attrs
            # self.data[key]['wigner'].attrs['dataType'] = 'wigner'
            # self.data[key]['wigner'].attrs['units'] = 'attosecond'

            # Using function
            self.data[key]['wigner'] = mfWignerDelay(self.data[key]['TX'], pType = pType)



    def phaseCorrection(self, keys = None, cPhase = None, lPhase = True, **kwargs):
        """
        Phase correction function - apply Coulomb (or other) phase corrections.

        Wraps :py:func:`epsproc.calc.phases.phaseCorrection()` for self.data[key] items.

        Results are set to self.data[key]['matE'] and self.data[key]['phaseCorr']  (this stores additional phase term only).

        Parameters
        ----------
        keys : str, int or list, optional, default = None
            If set, use only these datasets (keys).
            Otherwise run for all datasets in self.data.

        cPhase : Xarray, optional, default = None
            Coulomb phases to use.
            If None, compute these for input array using :py:func:`epsproc.calc.phases.coulombPhase()` function.
            **kwargs are also passed to :py:func:`epsproc.calc.phases.coulombPhase()`

        lPhase : bool, default = True
            Apply :math:`i^{-l}` term?
            Skip this term if false

        **kwargs : optional keyword args
            Passed to :py:func:`epsproc.calc.phases.coulombPhase()`

        TODO:

        - fix ugly .unstack('LM') included here, should use smart dim handling.
        - similarly have .swap_dims({'Eke':'Eh'}) to force E dims, should be smartly handled via dim checks. Can also just pass Elist here?

        """

        # Check/set keys
        keys = self._keysCheck(keys)

        # Loop over datasets & store output
        for key in keys:

            dataCorr, phaseCorr = phaseCorrection(self.data[key]['matE'].unstack('LM').swap_dims({'Eke':'Eh'}), cPhase = cPhase, lPhase = lPhase, **kwargs)

            # Propagate attrs
            dataCorr.attrs = self.data[key]['matE'].attrs
            phaseCorr.attrs = self.data[key]['matE'].attrs
            # TX.attrs['dataType'] = 'mfpad'

            self.data[key]['matE'] = dataCorr.stack({'LM':['l','m']}).swap_dims({'Eh':'Eke'})
            self.data[key]['phaseCorr'] = phaseCorr



    # **** Quick hack IO for dataarrays
