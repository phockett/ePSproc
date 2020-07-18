"""
ePSproc wfPlot

Functions for calculating and plotting wavefunctions from ePolyScat.

Makes use of pyvista/itkwidgets on the backend (see notes below).

15/07/20    v1, adapting from test notebook, pyVista_tests_070320_basics_working_150320.ipynb
            See also :py:func:`epsproc.vol.orbPlot.molOrbPlotter` for class development, use as template.

"""

import pyvista as pv
import numpy as np

from epsproc.util import orb3DCoordConv
from epsproc.IO import readOrb3D


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

    """

    def __init__(self, fileIn = None, fileBase = None, fType = '_Orb.dat'):
        """Init wfPlotter. Read file(s) & set-up grids."""

        # Set and read source file or directory
        self.fileIn = fileIn
        self.fileBase = fileBase
        self.fType = fType

        self.dataSet = readOrb3D(fileIn=fileIn, fileBase=fileBase, fType=fType)

        print(f"Read {len(self.dataSet)} wavefunction data files OK.")

        self.setGrid()

        print('*** Grids set OK')

        self.setData()

        print('*** Data set OK')

        print(self.vol)


    def setGrid(self, methodType = 2):
        """
        Set wavefunction grids from files.

        Basic wrapper for :py:func:`epsproc.conversion.orb3DCoordConv()`.

        """

        # Q: assume grid same for all files?
        # Split by files and/or symmetry and/or energy?
        # Separate pyVista objects for each, or can use for multiple?
        # Basically comes down to datastructures - list, dict, xarray...?

        #*** Option (1): Loop over data files & append PV objects
        if methodType == 1:
            for n, fileIn in enumerate(self.dataSet):
                X,Y,Z = orb3DCoordConv(fileIn)
                vol = pv.StructuredGrid(X, Z, Y)
                self.dataSet[n].extend(vol)

        elif methodType == 2:
            #*** Option (2): use just first file & set shared PV object
            X,Y,Z = orb3DCoordConv(self.dataSet[0])

        # Set pyVista object to hold orbital calculations.
        # Note ordering!
        self.vol = pv.StructuredGrid(X, Z, Y)


    def setData(self):
        """
        Sort self.dataSet items to pyVista object.

        TODO: useful naming for data + structure by E and Sym.

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
        for n, fileIn in enumerate(self.dataSet):
            self.vol.point_arrays[f'{str(n)}-Re']=fileIn[3][0].flatten(order="F")
            self.vol.point_arrays[f'{str(n)}-Im']=fileIn[3][1].flatten(order="F")
            self.vol.point_arrays[f'{str(n)}-Abs']=fileIn[3][2].flatten(order="F")



    def plotWf(self, wfN = None, pType = 'Abs', isoLevels = 6, isoValsAbs = None, isoValsPC = None, interactive = True, opacity = 0.5):
        """
        Plot wavefunction(s) with pyVista/ITK

        wfN : str or list of strs, optional, default = None
            Wavefn(s) to plot, by name in self.vol PV object.
            By default these are set to wavefunction numbers, corresponding to files read in.
            If not supplied, plot all available surfaces, i.e. items in self.vol.array_names

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
        if wfN is not None:
            if type(wfN) is str:
                wfN = list(wfN)

            # Append plot type
            wfN = [f'{item}-{pType}' for item in wfN]

        else:
            wfN = self.vol.array_names


        for item in wfN:
            # Set plotting by type
            if item.endswith(pType):
                # Set limits.
                limitVals = [self.vol[item].min(), self.vol[item].max(), np.abs(self.vol[item]).mean()]

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
                pl.add_mesh(self.vol.contour(isosurfaces = isoValsOrb, scalars = item), smooth_shading=True, opacity=opacity)  # Plot iso = 0.1

        # Render plot
        # In notebook tests this doesn't reneder unless called again?
        self.pl = pl
        self.pl.show()
