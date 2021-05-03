"""
ePSproc wfPlot

Functions for calculating and plotting wavefunctions from ePolyScat.

Makes use of pyvista/itkwidgets on the backend (see notes below).

15/07/20    v1, adapting from test notebook, pyVista_tests_070320_basics_working_150320.ipynb
            See also :py:func:`epsproc.vol.orbPlot.molOrbPlotter` for class development, use as template.

"""

# Plotting & numerics
import pyvista as pv
import numpy as np

# File & path handling
from pathlib import Path
import os
import tempfile
import inspect
import json

# For inline gif display in notebook - should add env check here first.
from IPython.display import Image

# Local functions
from epsproc.util import orb3DCoordConv
from epsproc.util.misc import fileListSort, timeStamp  # 31/07/20 just added and not propagated module changes as yet.
from epsproc.IO import readOrb3D, getFiles

from epsproc.vol.setOptions import setLocalOptions, readOptionsFile, writeOptionsFile  # Plot options IO




class wfPlotter():
    """
    Class to define & plot wavefunctions from ePS data.

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

    prefix : str, optional, default = None
        In case of dir/file parsing issue, a manually set file prefix can be passed.
        This may be neceassy if filename components are not split by '_', or only a single file is present.
        (See example notebook somewhere...?)

    mol : str, optional, default = None
        Will be (hopefully) found from filenames, but can be set here.

    fType : str, optional, default = '_Orb.dat'
        File ending for ePS output files containing wavefunction data.

    dataTypes : list, optional, default = ['Re', 'Im', 'Abs']
        Datatypes present in data file (will be looped over on read in).

    optionsFile : str or path object, optional, default = None
        Specify path to plot options file (.json) used by setPlotOptions() method.
        Default file should be `[modulePath]epsproc/vol/plotOptions.json`, this will be set if nothing is passed.

    tempdir : str or path object, optional, default = None
        Temporary directory for image file outputs.
        Default is system tempdir, as provided by the `tempfile` module (`tempfile.gettempdir()`).




    Returns
    -------
    wfPlotter class object, with wavefunction data + plotting methods.


    Notes
    ------
    On the backend, this uses:
    - :py:func:`epsproc.readOrb3D()` for file IO
    - pyvista/ITK for plotting.
    - Should improve filehandling with Path() here...? See orbPlot.py.

    May want to change method ordering and structure here - this currently imports data, then assigns (duplicates) to pyVista object. Not very efficient.

    TODO: add molecular structure... should be able to plot with orbPlot object...?

    """

    def __init__(self, fileIn = None, fileBase = None, prefix = None, mol = None,
                 fType = '_Orb.dat', dataTypes = ['Re', 'Im', 'Abs'],
                 optionsFile = None, tempdir = None):
        """Init wfPlotter. Read file(s) & set-up grids."""

        # Set and read source file or directory
        self.fileIn = fileIn
        self.fileBase = fileBase
        self.fType = fType
        self.dataTypes = dataTypes
        self.prefix = prefix   # Set for file sorting, can set manually to override default if required.
        self.mol = mol

        self.optionsFile = optionsFile  # Set to None, defaults will be set later.
        self.plotOptions = {}
        self.fDictSel = {}

        self.tempdir = tempdir

        self.getOrbFiles()

        self.setPlotOptions()  # Get plotter options from file.


    def getOrbFiles(self, verbose=True):
        """
        Call functionality from IO.getFiles to allow for more control (vs. readOrb3D, which will always get all files).

        Also sort outputs by (Sym, E) groups.

        """

        # Get full file list
        self.fList = getFiles(fileIn = self.fileIn, fileBase = self.fileBase, fType = self.fType)

        # Sort list by (Sym, E), and set values to output dictionaries
        self.sortOrbFiles(masterFlag = True)

        # Basic error check for number of E points per symmetry
        Elist = [len(self.fDictSym[item]['E']) for item in self.fDictSym]

        if len(np.unique(Elist)) > 1:
            print(f"*** Warning: inconsitent E points per symmetry, {Elist}")

        # Print results
        if verbose:
            print(f"Found molecule: {self.mol}")
            print(f"Found {len(self.fDictSym)} symmetries, {[self.fDictSym[item]['sym'] for item in self.fDictSym]}")
            print(f"Found {len(self.fDictSym[0]['E'])} energies, {self.fDictSym[0]['E']}")


    def sortOrbFiles(self, fList = None, masterFlag = False, groupByPrefix = True):
        """
        Sort input list of wavefunction files to dictionary by symmetry & energy.

        Parameters
        ----------

        fList : list, optional, default = None
            List of files to sort.
            If fList = None, use self.fList.

        masterFlag : bool, optional, default = True
            If masterFlag = True, set to class variables, otherwise return values to calling function.

        groupByPrefix : bool, optional, default = True
            Group sorted files by prefix?
            Will use self.prefix if set, otherwise determined from file names by :py:func:`epsproc.util.misc.fileListSort`.
            In latter case, self.prefix will be set if masterFlag = True

        """

        if fList is None:
            fList = self.fList

        # Sort files & group
        # if len(fList) > 1:
        #     # NOTE - this will return in E order if natsort lib is present, and usual filename convention.
        #     # May want to supplement with further sorting by E later?
        #     fList, groupedList, prefixStr = fileListSort(fList, groupByPrefix=groupByPrefix, prefixStr = self.prefix, verbose=True)
        #
        # # Ugly fix for single file passing.
        # # EDIT: Now fixed in fileListSort()
        # else:
        #     groupedList = fList
        #
        #     if self.prefix is None:
        #         prefixStr = Path(self.fileBase, Path(fList[0]).name.rsplit('_')[0]).as_posix()
        #     else:
        #         prefixStr = self.prefix
        #     # print(prefixStr)

        # NOTE - this will return in E order if natsort lib is present, and usual filename convention.
        # May want to supplement with further sorting by E later?

        fList, groupedList, prefixStr = fileListSort(fList, groupByPrefix=groupByPrefix, prefixStr = self.prefix, verbose=True)

        # Case for single file without preset self.prefix.
        if (self.prefix is None) and (prefixStr is None):
            prefixStr = Path(self.fileBase, Path(fList[0]).name.rsplit('_')[0]).as_posix()
        elif prefixStr is None:
            prefixStr = self.prefix



        # Parse list

        # Get molecule from file name - OK, but possibly issues with rolling into symmetry label here?
        # Check for local setting first, can then also use this to override
        # Assumes file name schema [mol][sym]_[XXeV]
        if hasattr(self, 'mol') and (self.mol is not None):
            mol = self.mol
        else:
            # mol = prefixStr.split('/')[-1][0:-1]
            mol = Path(prefixStr).parts[-1][0:-1]  # File system agnostic version

        # Loop over groups & extract details.
        fDictSym = {}
        fDictE = {}
        for n, item in enumerate(groupedList):

            if type(item) != list:  # Fix for case of single file item.
                item = [item]

            itemSorted, *_ = fileListSort(item, groupByPrefix=True, verbose=False)
            sym = item[0].replace(prefixStr,'S').split('_')[0]  # Take symmetry as first item here, assume 'S' is chopped!
            E = [float(fileItem.replace(prefixStr,'').split('_')[1].strip('eV')) for fileItem in itemSorted]

            # Store as dictionary by symmetry group
            fDictSym[n] = {'mol':mol,
                            'sym':sym,
                            'fList':itemSorted,
                            'E':E}

            # # List by E - now set below
            # for Ekey in E:
            #     if Ekey in fDictE.keys():
            #         fDictE[Ekey]['fList'].append(itemSorted)
            #     else:
            #         fDictE[Ekey]['fList'] = (itemSorted)

        # Dicionary sorted by E
        syms = [fDictSym[sym]['sym'] for sym in fDictSym]

        for n, Ekey in enumerate(E):
            fDictE[n] = {}
            fDictE[n]['E'] = Ekey
            fDictE[n]['fList'] = [fDictSym[sym]['fList'][n] for sym in fDictSym]
            fDictE[n]['mol'] = fDictSym[0]['mol']
            fDictE[n]['syms'] = syms

        # Set outputs to master (class) variables.
        if masterFlag:
            self.mol = mol
            self.E = E
            self.syms = syms
            self.prefix = prefixStr
            self.fList = fList
            self.fDictSym = fDictSym
            self.fDictE = fDictE

        else:
            return fList, fDictSym, fDictE


    def selectOrbFiles(self, fDictE = None, EList = None, SymList = None, verbose = True):
        """
        Subselect data files by Sym and E from sorted dict.

        Parameters
        ----------

        fDictE : dict, optional, default = None
            Pass dictionary for subselection.
            If not set, use self.fDictE.

        EList : list of ints, floats, optional, default = None
            Pass list of energies to select files to read in.
            If floats, these are assumed to be EKE values to match.
            If ints, these are used as indexes into the master EKE list.

        SymList : list of strs or ints, optional, default = None
            Pass list of symmetries to select for file IO.
            Either strings of symmetry names, or by index.

        verbose : bool, optional, default = True
            Print info on subselections if True.

        Returns
        -------
        self.fDictSel dictionary will be set to subselected items.


        Note
        -----

        Lists are parsed in order above, i.e. fDictE is set, then filtered by E and/or Sym.
        Unmatched inputs are ignored.


        TODO

        - Add a range option.
        - Nearest value matching.

        """


        # Subselect by Sym and E from dict.
        # Version with sym-sorted list - now use version with E sorting!
        # if SymList is not None:
        #     if type(SymList[0]) is str:
        #         fileDictSel = [fileDict[item] for item in fileDict if fileDict[item]['sym'] in SymList]
        #     else:
        #         # fileDictSel = [{item:fileDict[item]} for item in SymList]  # This will keep dict format, with extra [list] wrapper
        #         fileDictSel = [fileDict[item] for item in SymList]  # This will set to list, but for int indexing it's equivalent

        # Start with E-indexed dict
        # fileDictSel = self.fileDictE.copy()

        # Default to master list if not passed
        if fDictE is None:
            fDictE = self.fDictE.copy()


        if EList is not None:
            if type(EList[0]) is int:
                fileDictSel = {k: v for k, v in fDictE.items() if k in EList}  # For index list
            else:
                fileDictSel = {k: v for k, v in fDictE.items() if fDictE[k]['E'] in EList}
        else:
            # Default to full E-indexed dict
            fileDictSel = fDictE.copy()

        # For symmetries, subselect file and sym lists per key.
        # There's probably a neater way to do this, maybe with sets?
        if SymList is not None:
            if type(SymList[0]) is str:
                for key in fileDictSel.keys():
                    inds = [n for n,m in enumerate(fileDictSel[key]['syms']) if m in SymList]  # Get indexes
                    fileDictSel[key]['syms'] = [fileDictSel[key]['syms'][n] for n in inds]
                    fileDictSel[key]['fList'] = [fileDictSel[key]['fList'][n] for n in inds]
            else:
                inds = SymList
                for key in fileDictSel.keys():
                    fileDictSel[key]['syms'] = [fileDictSel[key]['syms'][n] for n in inds]
                    fileDictSel[key]['fList'] = [fileDictSel[key]['fList'][n] for n in inds]

        self.fDictSel = fileDictSel

        if verbose:
            print(f"Selected {len(fileDictSel)} energies, {[fileDictSel[item]['E'] for item in fileDictSel]}")
            print(f"Selected symmetries, {self.fDictE[list(self.fDictE.keys())[0]]['syms']}")


    def readOrbFiles(self, fList = None, EList = None, SymList = None, verbose = False):
        """
        Read wavefunction files.

        Parameters
        ----------
        flist : list of strs, Paths or ints, optional, default = None
            Pass fList as list of files or ints as indexes into self.fList index.
            Note this is the ungrouped file list (the original fList), and items are resorted, so may not match self.fDictE or self.fDictSym

        EList : list of ints, floats, optional, default = None
            Pass list of energies to select files to read in.
            If floats, these are assumed to be EKE values.
            If ints, these are used as indexes into the master EKE list.

        SymList : list of strs or ints, optional, default = None
            Pass list of symmetries to select for file IO.
            Either strings of symmetry names, or by index.

        Lists are parsed in order above, i.e. fList is set, then filtered by E and/or Sym.
        If nothing is set, read all files from self.fList.

        """

        # Set master file list.
        selFlag = False

        # Set cases for using self.fDictSel or passed args.
        # TODO: better logic here - switch to **kwargs?
        if not self.fDictSel:
            selFlag = True
        elif (fList is not None):
            selFlag = True
        elif (EList is not None):
            selFlag = True
        elif (SymList is not None):
            selFlag = True

        if selFlag:
            if fList is not None:
                if type(fList[0]) is str:
                    fileIn = fList
                elif type(fList[0]) is int:
                    # fileIn = self.fList[fList]  # Select items from existing list
                    fileIn = [self.fList[n] for n in fList]
                else:
                    print(f"File list item type {type(fList[0])} not supported.")
            else:
                fileIn = self.fList

            # Sort input list
            # print(fileIn)
            fileIn, fileDictSym, fileDictE = self.sortOrbFiles(fList=fileIn, masterFlag=False)

            # Subselect files
            self.selectOrbFiles(fDictE = fileDictE, EList = EList, SymList = SymList)

        # # Stack to list for file IO function.
        # # fileIn = [self.fDictSel[key]['fList'] for self.fDictSel.keys()]  # OK, but not flat list
        # fileIn = [item for key in self.fDictSel.keys() for item in self.fDictSel[key]['fList']]  # Flat list over keys.
        self.fListSel = [item for key in self.fDictSel.keys() for item in self.fDictSel[key]['fList']]  # Flat list over keys.

        # # Read files as set
        # print(f"Reading {len(fileIn)} wavefunction data files.")
        # # self.dataSet = readOrb3D(fileIn = fileIn, fileBase = self.fileBase, fType = self.fType)
        # dataSet = readOrb3D(fileIn = fileIn, fileBase = self.fileBase, fType = self.fType)
        # print(f"Read {len(self.dataSet)} wavefunction data files OK.")

        # Alternative method - set data in fDictSel to maintain traceability & labels.
        fTot = 0
        for key in self.fDictSel.keys():
            self.fDictSel[key]['dataSet'] = readOrb3D(fileIn = self.fDictSel[key]['fList'], fileBase = self.fileBase, fType = self.fType, verbose = verbose)
            fTot += len(self.fDictSel[key]['fList'])

        print(f"\nRead {fTot} wavefunction data files OK.")

        self.setGrid()

        print('*** Grids set OK')

        self.setData()

        print('*** Data set OK')

        # print(self.vol)



    def setGrid(self, methodType = 2):
        """
        Set wavefunction grids from files.

        Basic wrapper for :py:func:`epsproc.conversion.orb3DCoordConv()`.

        """

        # Q: assume grid same for all files?
        # Split by files and/or symmetry and/or energy?
        # Separate pyVista objects for each, or can use for multiple?
        # Basically comes down to datastructures - list, dict, xarray...?
        # 03/08/20 - now loops over items in fDictSel, and applies method to each key.

        for key in self.fDictSel:
            #*** Option (1): Loop over data files & append PV objects
            if methodType == 1:
                for n, fileIn in enumerate(self.fDictSel[key]['dataSet']):
                    X,Y,Z = orb3DCoordConv(fileIn)
                    vol = pv.StructuredGrid(X, Z, Y)
                    self.fDictSel[key]['dataSet'][n].extend(vol)

            elif methodType == 2:
                #*** Option (2): use just first file & set shared PV object
                X,Y,Z = orb3DCoordConv(self.fDictSel[key]['dataSet'][0])

                # Set pyVista object to hold orbital calculations.
                # Note ordering!
                self.fDictSel[key]['vol'] = pv.StructuredGrid(X, Z, Y)

        if methodType == 3:
            # Set one master array for *all* datasets
            X,Y,Z = orb3DCoordConv(self.fDictSel[self.fDictSel.keys()[0]]['dataSet'][0])
            self.vol = pv.StructuredGrid(X, Z, Y)


    def setData(self):
        """
        Sort self.dataSet items to pyVista object.

        TODO: useful naming for data + structure by E and Sym.
        EDIT: this is now in self.fDictSym, self.fDictE, self.fDictSel, and should be propagated here.

        NOTE - PyVista items are indexed sequentially, so may not match input dict keys. This may change in future.
        """

        # Currently set in init(), but may want to move here?
        # self.dataSet = readOrb3D(fileIn = self.fileIn, fileBase = self.fileBase, fType = self.fType)

        # Test method
        # wfPv.point_arrays['Re']=wfData[0][3][0].flatten(order="F")
        # wfPv.point_arrays['Im']=wfData[0][3][1].flatten(order="F")
        # wfPv.point_arrays['Abs']=wfData[0][3][2].flatten(order="F")

        # Updated for this code
        # NOPE - complex values not supported. Instead, set at plot time?
        # for n, fileIn in enumerate(self.dataSet):
        #     self.vol.point_arrays[str(n)] = (fileIn[3][0] + 1j*fileIn[3][1]).flatten(order="F")

        # In this case set all data.
        # May be better to set only plot data in plotting fn?
        # for n, fileIn in enumerate(self.dataSet):
        #     self.vol.point_arrays[f'{str(n)}-Re']=fileIn[3][0].flatten(order="F")
        #     self.vol.point_arrays[f'{str(n)}-Im']=fileIn[3][1].flatten(order="F")
        #     self.vol.point_arrays[f'{str(n)}-Abs']=fileIn[3][2].flatten(order="F")

        # Version for dict structure

        # Set Re, Im and Abs arrays
        arrayInputs = self.dataTypes

        # Set global outputs
        globalLimits = {}
        for item in arrayInputs:
            globalLimits[item] = []

        # Loop over files, and set to vol
        for key in self.fDictSel:
            for n, fileIn in enumerate(self.fDictSel[key]['dataSet']):

                # Loop over input data and assign to Pyvista object.
                for m, item in enumerate(arrayInputs):
                    arrayName = f"{key}_E{self.fDictSel[key]['E']}_{self.fDictSel[key]['syms'][n]}-{item}"
                    self.fDictSel[key]['vol'].point_arrays[arrayName]=fileIn[3][m].flatten(order="F")

                    # Set also a dict of item properties per key, for later reference
                    self.fDictSel[key][item] = {'limits':[self.fDictSel[key]['vol'].point_arrays[arrayName].min(),
                                                            self.fDictSel[key]['vol'].point_arrays[arrayName].max(),
                                                            np.abs(self.fDictSel[key]['vol'].point_arrays[arrayName]).mean()]}


                    globalLimits[item].append(self.fDictSel[key][item]['limits'])

                # Basic version without inner loop
                # self.fDictSel[key]['vol'].point_arrays[f"{key}_E{self.fDictSel[key]['E']}_{self.fDictSel[key]['syms'][n]}-Re"]=fileIn[3][0].flatten(order="F")
                # self.fDictSel[key]['vol'].point_arrays[f"{key}_E{self.fDictSel[key]['E']}_{self.fDictSel[key]['syms'][n]}-Im"]=fileIn[3][1].flatten(order="F")
                # self.fDictSel[key]['vol'].point_arrays[f"{key}_E{self.fDictSel[key]['E']}_{self.fDictSel[key]['syms'][n]}-Abs"]=fileIn[3][2].flatten(order="F")

        # Set global limits
        self.globalLimits = {}
        for item in arrayInputs:
            self.globalLimits[item] = np.asarray(globalLimits[item])


    def head(self):
        """Return first item in self.fDictSel[key]"""

        itemKey = list(self.fDictSel.keys())[0]
        print(f"Item key: {itemKey}")

        return self.fDictSel[list(self.fDictSel.keys())[0]]

    def tail(self):
        """Return last item in self.fDictSel[key]"""
        itemKey = list(self.fDictSel.keys())[-1]
        print(f"Item key: {itemKey}")

        return self.fDictSel[list(self.fDictSel.keys())[-1]]

    def writeOptions(self):
        """Wrapper for file writer"""

        if self.optionsFile is not None:
            writeOptionsFile(self.optionsFile, self.plotOptions)

        else:
            print("No file set for plot option file writer.")


    def setPlotOptions(self, optionsFile = None, verbose = False, **kwargs):   #, recursive = True):
        """
        List of default plot options. Set here, change later if required.

        Default file should be `[modulePath]epsproc/vol/plotOptions.json`, this will be set if nothing is passed.

        TODO: testing to see if this is robust. Otherwise may want to include defaults directly in source here.
        Setting options in nested dicts from **kwargs needs some work.
        See: http://127.0.0.1:8888/notebooks/python/python_basics/Dictionary_methods_Benedict_120820.ipynb

        """

        # Default values
        if optionsFile is None:
            if self.optionsFile is not None:
                optionsFile = self.optionsFile

            else:
                # Set path based on file location - may not be robust?
                # From https://stackoverflow.com/a/12154601
                optionsFile = Path((os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))),'plotOptions.json')
                self.optionsFile = optionsFile

        # Set path wrapper in case str was passed.
        # optionsFile = Path(optionsFile)

        # Call file read function
        owFlag = True
        if self.plotOptions and not kwargs:
            owFlag = (True if input("Overwrite local plot options with settings from file (y/n)? ") == 'y' else False)

        if self.plotOptions and kwargs:
            if verbose:
                print('Adding passed args to plotOptions')

            print('Arb option setting not yet supported, but can be manually assigned to self.plotOptions dictionary.')
            # for key in kwargs.keys():

            owFlag = False


        if owFlag:
            if verbose:
                print(f'Reading options from file {optionsFile}')

            self.plotOptions = readOptionsFile(optionsFile = optionsFile, verbose = verbose)

        #
        # if optionsFile.is_file():
        #     with open(optionsFile) as json_file:
        #         optionsFileJSON = json.load(json_file)
        #
        #     print(f"\n*** Read existing plot options from file {optionsFile} OK.")
        #     if verbose:
        #         print(json.dumps(optionsFileJSON, sort_keys=False, indent=4))
        #
        #     self.plotOptions = optionsFileJSON
        #
        # else:
        #     print(f"\n*** Plot options file {optionsFile} not found, using defaults.")
        #
        #     self.plotOptions = setLocalOptions()

            # if recursive = True:
            #     self.plotOptions = None  # TODO - set defaults here.
            # else:
            #     pass  # Rewrite options file here


        # # Default values
        # if (optionsFile is None) or (self.plotOptions is None):
        #     # Set path based on file location - may not be robust?
        #     # From https://stackoverflow.com/a/12154601
        #     optionsFile = Path((os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))),'plotOptions.json')
        #     # Call method recursively with new optionsFile.
        #     # Set recursive = False to rewrite the file if it doesn't exist this time!
        #     self.setPlotOptions(optionsFile = optionsFile, recursive = False)


    def plotWf(self, pType = None, opacity = None, verbose = True):
        """
        Plot wavefunction(s) with pyVista/ITK.

        ALREADY BECOMING SPAGHETTI... must do better here.

        07/08/20    v3, use this as main wrapper & to set logic, then call sub-methods for various plotting cases.
                        Options to a dict, this should fix issues with some methods only available for some plot types.

        TODO
        - kwargs pass to setLocalOptions and/or pyvista?
        - Logic/hierarchy for plot type and options. Setting plotter-based options may help?
        - Orbit animations, combined with subplots?

        """


        # Set plotter based on options
        if self.plotOptions['global']['inline']:
            if self.plotOptions['global']['subplot']:


                if self.plotOptions['global']['animate']:
                    # pl = pv.Plotter(shape = (2,2))  # Set for 2x2 grid with slice views
                    pl = pv.Plotter(shape = (1,2)) # Set for 1x2 grid with slice views

                # Set for subplot per dataset
                else:
                    pN = len(self.fDictSel)
                    rP = np.rint(np.sqrt(pN)).astype(np.int32)
                    cP = np.ceil(pN/rP).astype(np.int32)
                    sPlots=(rP,cP)
                    currPlot = [0, 0]
                    pl = pv.Plotter(shape = sPlots)

                # Temporary hack - reset interactive here to avoid plotter related issues.
                self.plotOptions['global']['interactive'] = False

            elif self.plotOptions['global']['interactive'] and (not self.plotOptions['global']['subplot']):
                pl = pv.PlotterITK()

            else:
                pl = pv.Plotter()
                # For animation case may want to add screen size and off_screen?
                # https://github.com/pyvista/pyvista-support/issues/81#issuecomment-563493574

        else:
            pl = pv.BackgroundPlotter()


        self.pl = pl

        # Set up animation stuff
        if self.plotOptions['global']['animate']:

            fN = 0  # Frame counter
            fString = f'wfAnimation_{self.mol}_{timeStamp()}' # File name schema

            # Set temp dir - set manually to allow for persistence
            if self.tempdir is None:
                self.tempdir = tempfile.gettempdir()  # Set to system default if not supplied

            tempdir = Path(self.tempdir, fString)
            tempdir.mkdir()

            print('Frame images output to: ', tempdir)
            # print('TODO: animate these.')
            frames = []

            self.gifFile = tempdir.joinpath(f'{format(fString)}.gif')
            self.pl.open_gif(self.gifFile.as_posix())
            print('GIF animation file: ', self.gifFile)


        if verbose:
            print(f"Set plotter to {pl.__class__.__name__}")

        # Parse additional passed args
        if opacity is None:
            opacity = self.plotOptions['global']['opacity']
        else:
            self.plotOptions['global']['opacity'] = opacity

        # Set iso vals for plot, per dataset or globally
        self.setIsoVals()

        # Loop over data d add meshes
        # May want to move to a sub-function to help with animation cases.
        if pType is None:
            pType = self.plotOptions['global']['pType']
        else:
            self.plotOptions['global']['pType'] = pType  # Update options if set.


        for key in self.fDictSel:
            wfN = self.fDictSel[key]['vol'].array_names
            vol = self.fDictSel[key]['vol']   # Set pointer here for brevity

            for item in wfN:
                if item.endswith(pType):

                    # Calculate contours with selected scalars
                    contourItem = vol.contour(isosurfaces = self.fDictSel[key][pType]['isoValsOrb'],scalars = item)

                    # TODO - add subplot options - now set for slices below
                    # if self.plotOptions['global']['subplot']:
                    #     pass

                    # For animation save frame then clear
                    if self.plotOptions['global']['animate']:

                        # 29/04/21 - added multiple slice planes subplotting option, based on https://docs.pyvista.org/examples/02-plot/ortho-slices.html#sphx-glr-examples-02-plot-ortho-slices-py
                        if self.plotOptions['global']['subplot']:
                            # slices = contourItem.slice_orthogonal(x=0, y=0, z=0)  # Slice from contours
                            slices = vol.slice_orthogonal(x=0, y=0, z=0)  # Slice from full volume


                            # ADditional options
                            #dargs = dict(cmap='gist_ncar_r')
                            dargs = {}

                            # Plot 2x2 grid
                            # p = pv.Plotter(shape=(2,2))
                            # XYZ - show 3D scene first
                            # self.pl.subplot(1,1)
                            self.pl.subplot(0,0)
                            # p.add_mesh(slices, **dargs)  # As slices
                            self.pl.add_mesh(contourItem, smooth_shading=True, opacity=opacity, clim=self.fDictSel[key][pType]['cmapLimitVals'])  # As vol
                            self.pl.show_grid()
                            # p.camera_position = cpos
                            self.pl.view_isometric()
                            self.pl.add_text("{: >8.2f} eV".format(self.fDictSel[key]['E']), position='upper_right')

                            # # XY
                            # self.pl.subplot(0,0)
                            # self.pl.add_mesh(slices, **dargs)
                            # self.pl.show_grid()
                            # self.pl.camera_position = 'xy'
                            # self.pl.enable_parallel_projection()
                            # # ZY
                            # self.pl.subplot(0,1)
                            # self.pl.add_mesh(slices, **dargs)
                            # self.pl.show_grid()
                            # self.pl.camera_position = 'zy'
                            # self.pl.enable_parallel_projection()
                            # # XZ
                            # self.pl.subplot(1,0)
                            # self.pl.add_mesh(slices, **dargs)
                            # self.pl.show_grid()
                            # self.pl.camera_position = 'xz'
                            # self.pl.enable_parallel_projection()

                            self.pl.subplot(0,1)
                            self.pl.add_mesh(slices, **dargs)
                            self.pl.show_grid()
                            self.pl.camera_position = 'yz'
                            self.pl.enable_parallel_projection()

                        else:

                            self.pl.add_mesh(contourItem, smooth_shading=True, opacity=opacity, clim=self.fDictSel[key][pType]['cmapLimitVals'])  # Plot iso = 0.1
                            # cpos = self.pl.render(cpos = cpos)  # TODO - camera position setting & matching, orbits & pans.
                            # self.pl.add_axes()  # Throws error ERROR:root:The interactor must be set prior to enabling/disabling widget
                            self.pl.view_isometric()  # Add this explicitly, otherwise get a flat (x,y) view.
                            # self.pl.add_text(f"{self.fDictSel[key]['E']} eV", position='lower_left')
                            self.pl.add_text("{: >8.2f} eV".format(self.fDictSel[key]['E']), position='upper_right')

                        # Following example code... but not quite correct here.
                        # For now just use render() option.
                        # if fN == 0:
                        #     self.pl.show(auto_close=False)
                        # else:
                        #     self.pl.render()

                        self.pl.render()

                        frames.append(self.pl.screenshot(filename = tempdir.joinpath("f{:04}_{}.png".format(fN,fString)).as_posix(), return_img = True))
                        self.pl.write_frame()

                        if self.plotOptions['global']['subplot']:
                            # self.pl.subplot(1,1)  # Clear iso plot specifically.
                            self.pl.subplot(0,0)
                            self.pl.clear()
                        else:
                            self.pl.clear()
                        fN += 1

                    # Otherwise plot with/without subplotting
                    else:
                        if self.plotOptions['global']['subplot']:
                            # Set
                            self.pl.subplot(currPlot[0], currPlot[1])
                            # Increment
                            currPlot[1] += 1
                            if currPlot[1] > (sPlots[1]-1):  # Wrap rows - probably a neater way to do this...!
                                currPlot[1] = 0
                                currPlot[0] += 1

                            # if np.base_repr(currPlot[1], base=sPlot[1]) == 10:  # Alternative with base set.


                        # Add contours for currently selected scalars (orbital)
                        if self.plotOptions['global']['interactive'] and hasattr(self.pl,'add_mesh_clip_plane'):
                            self.pl.add_mesh_clip_plane(contourItem, smooth_shading=True, opacity=opacity)
                        else:
                            self.pl.add_mesh(contourItem, smooth_shading=True, opacity=opacity)  # Plot iso = 0.1


        # Additional options (for some plotters)
        if hasattr(self.pl, 'add_axes'):
            self.pl.add_axes()

        # TODO:
        # Clipping plane?
        # plotter.add_mesh_clip_plane
        # Subplots option
        # ANIMATION

        # Write gif
        if self.plotOptions['global']['animate']:
            self.animationPath = [tempdir, fString]
            self.frames = frames

            self.pl.close()

            print('Animation output complete.')
        #     from PIL import Image, ImageDraw
        #     Image.save('out.gif', save_all=True, append_images=frames)

            if self.plotOptions['global']['inline']:
                # Display gif (in notebook)
                # From https://stackoverflow.com/a/59988258
                # Image(self.gifFile)

                # More sophisticated version
                return self.showGif()

        elif self.plotOptions['global']['inline']:
            # Setting the return value here shows the plot in a notebook.
            # Otherwise, need to call self.pl.show() in notebook after running method.
            return self.pl.show()

    # More sophisticated version from
    # https://github.com/ipython/ipython/issues/10045#issuecomment-642640541
    # This ensures HTML export from Jupyter notebooks displays OK
    def showGif(self):
        """
        Show GIF in Jupyter notebook usng IPython.display

        HTML export safe version, from https://github.com/ipython/ipython/issues/10045#issuecomment-642640541
        Thanks to @maedoc Marmaduke Woodman, https://github.com/maedoc

        """

        import base64
        from IPython import display

        with open(self.gifFile, 'rb') as fd:
            b64 = base64.b64encode(fd.read()).decode('ascii')

        return display.HTML(f'<img src="data:image/gif;base64,{b64}" />')


    def setIsoVals(self):
        """
        Set isosurf properties per item, or use global values.
        """

        # Loop over data and set isosurf properties
        for dataType in self.dataTypes:
            # Set global limits
            limitValsGlobal = self.globalLimits[dataType].max(axis=0)

            for key in self.fDictSel:

                # Select global or per-array limits
                if self.plotOptions['global']['isoValsGlobal']:
                    limitVals = limitValsGlobal
                else:
                    limitVals = self.fDictSel[key][dataType]['limits']

                # Set default case
                isoValsOrb = np.linspace(-limitVals[2], limitVals[2], self.plotOptions['global']['isoLevels'])  # Again set from mean.

                # Override if alternative limits set (logic may be dodgy here)
                if self.plotOptions['global']['isoValsAbs'] is not None:
                    isoValsOrb = np.array(self.plotOptions['global']['isoValsAbs'])

                if self.plotOptions['global']['isoValsPC'] is not None:
                    isoValsPC = np.array(self.plotOptions['global']['isoValsPC'])
                    isoValsOrb = np.r_[isoValsPC, -isoValsPC] * limitVals[2]  # Set PC vals from mean?

                self.fDictSel[key][dataType]['isoValsOrb'] = isoValsOrb
                self.fDictSel[key][dataType]['cmapLimitVals'] = [isoValsOrb[0], isoValsOrb[-1]] # Set to use for global cmapping


    def plotWfV2(self, wfN = None, pType = 'Abs', isoLevels = 6, isoValsAbs = None, isoValsPC = None,
                interactive = True, opacity = 0.5, animate = False):
        """
        Plot wavefunction(s) with pyVista/ITK

        # CURRENTLY NOT FUNCTIONAL with dictionary version.
        # SUBSELECT files first instead.
        # wfN : str or list of strs, optional, default = None
        #     Wavefn(s) to plot, by key (item number) from self.fDiictSel[key].vol PV object.
        #     By default these are set to wavefunction numbers, corresponding to files read in.
        #     If not supplied, plot all available surfaces, i.e. items in self.vol.array_names

        isoLevels : int, optional, default = 6
            Number of isosurfs to compute & render.

        isoValsAbs : list or array, default = None
            Isovals to use (absolute values).
            If both isoValsAbs and isoValsPC are passed, only PC values will be used.

        isoValsPC : list or array, default = None
            Isovals to use, percentage values. Plot values are set by +/-isoValsPC * np.abs(self.vol[item]).mean().
            If both isoValsAbs and isoValsPC are passed, only PC values will be used.

        interactive : bool, optional, default = True
            Plot inteactively using itkwidgets if true.
            Otherwise use pyVista.Plotter(), which produces static output.

        opacity : float, optional, default = 0.5
            Set isosurface opacity.

        animate : bool, optional, default = False
            Generate animation from dataset, as per https://docs.pyvista.org/examples/02-plot/gif.html
            Note this overrides the "interactive" setting above, and will only return the final frame to the self.pl object.
            UPDATE 06/08/20: plotter.update() for pv.Plotter() object doesn't seem to work for isosurfs, employing alternative method from https://docs.pyvista.org/plotting/plotting.html#plot-time-series-data
            This means that animation runs will use BackgroundPlotter() plotter, and open in a separate window.

        Notes
        -----
        07/08/20    v3, unthreading looped plotting in favour of two-method soloution - currently set as new method.
        03/08/20    v2, testing animation plus additional plotting options.
        18/07/20    v1, adapted from molOrbPlotter.plotOrb function.

        To do
        -----
        - Plot type (re/im/abs) to add. Set to abs for testing.
        - Fix naming & colourmapping for multiple objects in ITK plotter.
        - Consolidate atoms > molecular geometry pyVista object.
        - orbN as list of ints or np.array? Integrate with calcOrb()?
        - Opacity mapping for isosurfs? For pv.Plotter() can do this via transfer fns, e.g. 'linear', https://docs.pyvista.org/examples/02-plot/opacity.html#sphx-glr-examples-02-plot-opacity-py

        """

        # If  animation set, override interactive setting
        fN = 0  # Frame counter
        if animate:
            # Imports for testing - to be moved.
            # NOTE - these are used in https://docs.pyvista.org/plotting/plotting.html#plot-time-series-data
            # But not required for non-interactive plot-to-file case.
            # from threading import Thread
            # import time

            interactive = False
            fileOut = f"wfAnimation_{self.mol}_{timeStamp()}.gif"
            print(f"Animating frames, output file: {fileOut}")

            # Follow example at https://docs.pyvista.org/plotting/plotting.html#plot-time-series-data
            # Use BG plotter, and threaded output.
            pl = pv.BackgroundPlotter()

        else:
            # Set ITK widgets or pv.Plotter
            # May also want to add other methods here, e.g. pv.BackgroundPlotter()  # Runs plotter in separate window process.
            if interactive:
                pl = pv.PlotterITK()
            else:
                pl = pv.Plotter()



        # Set meshes to plot, either passed or all items in self.vol
        # if wfN is not None:
        #     if type(wfN) is str:
        #         wfN = list(wfN)
        #
        #     # Append plot type
        #     wfN = [f'{item}-{pType}' for item in wfN]
        #
        # else:
        #     wfN = self.vol.array_names
        # 03/08/20 - dict version, now set below. Need to add selection logic back in here.

        for key in self.fDictSel:

            wfN = self.fDictSel[key]['vol'].array_names
            vol = self.fDictSel[key]['vol']   # Set pointer here for brevity

            for item in wfN:
                # Set plotting by type
                if item.endswith(pType):
                    # Set limits.
                    # limitVals = [vol[item].min(), vol[item].max(), np.abs(vol[item]).mean()]

                    # Set default case
                    isoValsOrb = np.linspace(-limitVals[2], limitVals[2], isoLevels)  # Again set from mean.

                    # Override if alternative limits set (logic may be dodgy here)
                    if isoValsAbs is not None:
                        isoValsOrb = np.array(isoValsAbs)

                    if isoValsPC is not None:
                        isoValsPC = np.array(isoValsPC)
                        isoValsOrb = np.r_[isoValsPC, -isoValsPC] * limitVals[2]  # Set PC vals from mean?

                    # print(isoValsOrb)


                    # Additional settings for animation
                    if animate:
                        # Add label
                        pl.add_text(f"{self.fDictSel[key]['E']} eV", position='lower_left')

                        if fN == 0:
                            # Plot first item
                            pl.add_mesh(vol.contour(isosurfaces = isoValsOrb, scalars = item), smooth_shading=True, opacity=opacity)  # Plot iso = 0.1

                            # setup camera and close
                            pl.show(auto_close=False)

                            # Init file write
                            pl.open_gif(fileOut)

                        else:
                            # pl.update_coordinates(pts)
                            # pl.update_scalars(vol.contour(isosurfaces = isoValsOrb, scalars = item), smooth_shading=True, opacity=opacity)
                            # pl.update_scalars(vol.contour(isosurfaces = isoValsOrb, scalars = item))  # Keywords not accepted for update fn?
                            pl.update(vol.contour(isosurfaces = isoValsOrb, scalars = item))

                        pl.write_frame()
                        fN += 1

                    else:
                        # Add contours for currently selected scalars (orbital)
                        pl.add_mesh(vol.contour(isosurfaces = isoValsOrb, scalars = item), smooth_shading=True, opacity=opacity)  # Plot iso = 0.1



        # Render plot
        # In notebook tests this doesn't reneder unless called again? (For ITK widgets case, but not for native pv.Plotter())
        if animate:
            pl.close()

        else:
            self.pl = pl
            self.pl.show()
