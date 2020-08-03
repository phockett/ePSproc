"""
ePSproc wfPlot

Functions for calculating and plotting wavefunctions from ePolyScat.

Makes use of pyvista/itkwidgets on the backend (see notes below).

15/07/20    v1, adapting from test notebook, pyVista_tests_070320_basics_working_150320.ipynb
            See also :py:func:`epsproc.vol.orbPlot.molOrbPlotter` for class development, use as template.

"""

import pyvista as pv
import numpy as np
from pathlib import Path

from epsproc.util import orb3DCoordConv
from epsproc.util.misc import fileListSort  # 31/07/20 just added and not propagated module changes as yet.
from epsproc.IO import readOrb3D, getFiles




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

    fType : str, optional
        File ending for ePS output files, default '_Orb.dat'


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

    def __init__(self, fileIn = None, fileBase = None, prefix = None, mol = None, fType = '_Orb.dat'):
        """Init wfPlotter. Read file(s) & set-up grids."""

        # Set and read source file or directory
        self.fileIn = fileIn
        self.fileBase = fileBase
        self.fType = fType
        self.prefix = prefix   # Set for file sorting, can set manually to override default if required.
        self.mol = mol

        self.getOrbFiles()


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
        if len(fList) > 1:
            # NOTE - this will return in E order if natsort lib is present, and usual filename convention.
            # May want to supplement with further sorting by E later?
            fList, groupedList, prefixStr = fileListSort(fList, groupByPrefix=groupByPrefix, prefixStr = self.prefix, verbose=True)

        # Ugly fix for single file passing.
        # EDIT: Now fixed in fileListSort()
        else:
            groupedList = fList

            if self.prefix is None:
                prefixStr = Path(self.fileBase, Path(fList[0]).name.rsplit('_')[0]).as_posix()
            else:
                prefixStr = self.prefix
            # print(prefixStr)


        # Parse list

        # Get molecule from file name - OK, but possibly issues with rolling into symmetry label here?
        # Check for local setting first, can then also use this to override
        # Assumes file name schema [mol][sym]_[XXeV]
        if hasattr(self, 'mol') and (self.mol is not None):
            mol = self.mol
        else:
            mol = prefixStr.split('/')[-1][0:-1]

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


    def readOrbFiles(self, fList = None, EList = None, SymList = None):
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
            self.fDictSel[key]['dataSet'] = readOrb3D(fileIn = self.fDictSel[key]['fList'], fileBase = self.fileBase, fType = self.fType)
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
        for key in self.fDictSel:
            for n, fileIn in enumerate(self.fDictSel[key]['dataSet']):
                self.fDictSel[key]['vol'].point_arrays[f"{key}_E{self.fDictSel[key]['E']}_{self.fDictSel[key]['syms'][n]}-Re"]=fileIn[3][0].flatten(order="F")
                self.fDictSel[key]['vol'].point_arrays[f"{key}_E{self.fDictSel[key]['E']}_{self.fDictSel[key]['syms'][n]}-Im"]=fileIn[3][1].flatten(order="F")
                self.fDictSel[key]['vol'].point_arrays[f"{key}_E{self.fDictSel[key]['E']}_{self.fDictSel[key]['syms'][n]}-Abs"]=fileIn[3][2].flatten(order="F")


    def plotWf(self, wfN = None, pType = 'Abs', isoLevels = 6, isoValsAbs = None, isoValsPC = None, interactive = True, opacity = 0.5):
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

        Notes
        -----
        18/07/20    v1, adapted from molOrbPlotter.plotOrb function.

        To do
        -----
        - Plot type (re/im/abs) to add. Set to abs for testing.
        - Fix naming & colourmapping for multiple objects in ITK plotter.
        - Consolidate atoms > molecular geometry pyVista object.
        - orbN as list of ints or np.array? Integrate with calcOrb()?
        - Opacity mapping for isosurfs? For pv.Plotter() can do this via transfer fns, e.g. 'linear', https://docs.pyvista.org/examples/02-plot/opacity.html#sphx-glr-examples-02-plot-opacity-py

        """

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
            vol = self.fDictSel[key]['vol']   # Set pointed here for brevity

            for item in wfN:
                # Set plotting by type
                if item.endswith(pType):
                    # Set limits.
                    limitVals = [vol[item].min(), vol[item].max(), np.abs(vol[item]).mean()]

                    # Set default case
                    isoValsOrb = np.linspace(-limitVals[2], limitVals[2], isoLevels)  # Again set from mean.

                    # Override if alternative limits set (logic may be dodgy here)
                    if isoValsAbs is not None:
                        isoValsOrb = np.array(isoValsAbs)

                    if isoValsPC is not None:
                        isoValsPC = np.array(isoValsPC)
                        isoValsOrb = np.r_[isoValsPC, -isoValsPC] * limitVals[2]  # Set PC vals from mean?

                    # print(isoValsOrb)

                    # Add contours for currently selected scalars (orbital)
                    pl.add_mesh(vol.contour(isosurfaces = isoValsOrb, scalars = item), smooth_shading=True, opacity=opacity)  # Plot iso = 0.1

        # Render plot
        # In notebook tests this doesn't reneder unless called again? (For ITK widgets case, but not for native pv.Plotter())
        self.pl = pl
        self.pl.show()
