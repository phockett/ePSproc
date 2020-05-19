"""
ePSproc orbPlot

Functions for calculating and plotting orbitals from quantum chemistry results.

Makes use of pyvista/itkwidgets, cclib and chemview/pyquante on the backend (see notes below).

13/05/20    v1, adapting from test notebook, cclib_chemlab_ITK_pyVista_dev_150420.ipynb

"""

# Required imports
import numpy as np
import sys
from pathlib import Path
import pyvista as pv
import cclib

# Set chemlab fn. import
def importChemLabQC(chemPath = None):
    """
    Import chemlab. This is set as optional, and a local path can be defiend, since it is deprecated and may not install.

    The chemlab.qc.wavefunction() method essentially repackages pyQuante functionality, for 3D volume calculation for orbitals defined by Gaussian basis sets.

    Code: https://github.com/chemlab/chemlab/blob/master/chemlab/qc/wavefunction.py

    Method: https://gabrielelanaro.github.io/blog/2015/01/14/visualizing-molecular-orbitals-ipython-notebook.html

    To do: add more error parsing and logging, see, e.g., https://airbrake.io/blog/python-exception-handling/importerror-and-modulenotfounderror

    """

    qcImport = False

    try:
        # import chemlab
        from chemlab.qc import wavefunction
        qcImport = True
    except ImportError:
        pass

    if chemPath is not None:
        try:
            # Load dev scripts - sidesteps issues with chemlab installation
            sys.path.append(chemPath)
            # from chemlab.qc import wavefunction  # Can't seem to get this version to work - likely due to init files?
            from qc import wavefunction  # This is OK with Base or subdir set
            qcImport = True
        except ImportError:
            pass

    if not qcImport:
        print('Chemlab orbital functionality not found, orbital calculations not possible.')
        return None
    else:
        print(f'Import OK: Chemlab module {wavefunction}')
        return wavefunction


class molOrbPlotter():
    """
    Class to define & calculate electronic structure properties & orbitals.

    Parameters
    ----------
    chemPath : str, optional
        Set import path for chemlib. If not set, only standard python dir tree will be used.

    fileIn : str, Path or jobInfo dict
        Set file to use, either directly via str/Path, or use file from jobInfo dictionary (as returned by :py:func:`epsproc.headerFileParse`)

    localPath : str, optional
        Set local path to use for electronic structure file, instead of original path from jobInfo.

    fileType : str, optional
        Set to override electronic structure file type by suffix.
        Usually this will be a .molden file for ePS, with a corresponding original file (e.g. Gamess .log).
        CClib supports most formats, although doesn't read Molden. See https://cclib.github.io/data.html

    Returns
    -------
    molOrbPlotter class object, with molecular data + plotting methods.

    To do
    -----
    - File serach logic? Set default case to call :py:func:`epsproc.headerFileParse`? Other locations?

    Electronic structure
    --------------------
    Read files and process with cclib.
    https://cclib.github.io/how_to_parse.html

    Orbitals
    --------
    Method from chemview.viewer.MoleculeViewer.add_isosurfaces
    The chemlab.qc.wavefunction() method essentially repackages pyQuante functionality, for 3D volume calculation for orbitals defined by Gaussian basis sets.

    """

    def __init__(self, chemPath = None, fileIn = None, localPath = None, fileType = None, geomInd = -1):
        """Init molOrbPlotter. Read electronic structure file & set-up grids."""

        self.wavefunction = importChemLabQC(chemPath = chemPath)

        if type(fileIn) is dict:
            # Get file name from jobInfo
            self.fileIn = Path(fileIn['Convert'][0].split()[1].strip("'"))

            # Override path if local path specified
            if localPath is not None:
                self.fileIn = Path(localPath, self.fileIn.name)

        else:
            self.fileIn = Path(fileIn)

        if fileType is not None:
            self.fileIn.suffix

        if self.fileIn.is_file():
            print(f'Found electronic structure file: {self.fileIn}')
            self.data = cclib.io.ccread(self.fileIn.as_posix())
            # print(self.data.metadata)
            print("Read %i atoms and %i MOs" % (self.data.natom, self.data.nmo))

            # Set which molecular geometry to use.
            self.geomInd = geomInd

            # Set dict for storing calc fns.
            self.f = {}

            # Set grid
            self.setGrid()

            # Set dictionary to hold volume calcs.
            # self.vol = {}

            # Set pyVista object to hold orbital calculations. NO SET IN self.setGrid()
            # Note ordering!
            # self.pvObj = pv.StructuredGrid(self.X, self.Z, self.Y) # OK

            print('*** Grids set OK')
            print(self.vol)

        else:
            print('***Missing electronic structure file: {self.filleIn}')




    def setGrid(self, extent = None, resolution = 50, padGrid = 1.0, sqGrid = True, geomInd = None):
        """Setup grids for orbital calculations.

        Basic 3D grid setup, based on passed or atomic coords.
        Currently limited to square grid for compatibility with ITK/VTK plotting - but should be able to make this more flexible.
        Final grid set as PyVista object, self.vol

        Parameters
        ----------
        extent : list or np.array, [min, max], optional
            Define limits of grid (in working units, Angs or Bohr).
            If not specified will be defined from atomic coords (min and max values in any dimension).

        resolution : int, optional, default = 50
            Number of points per dimension, final grid size will be resolution^3

        padGrid : float, optional, default = 1.0
            Padding to add when using atomic coords to define grid (in working units angs/bohr).
            Currently single valued, same for all dims.

        sqGrid : bool, optional, default = True
            Force square gridding? Currently only True is supported.

        geomInd : int, optional, default = -1
            Set to select a specific geometry from the electronic structure file (if present).
            Defaults to the final geometry listed.
            Overrides self.geomInd if set.

        Notes
        ------
        - Method originally adapted from chemview.viewer.MoleculeViewer.add_isosurfaces
        - Q: do units matter here...?

        """

        #*** Check global vals
        if geomInd is not None:
            self.geomInd = geomInd

        #*** Check and set grid limits
        if extent is None:
            # This will generate a [2x3] array of (min,max) values for (x,y,z) dims.
            extent = np.array([self.data.atomcoords[self.geomInd].min(axis=0), self.data.atomcoords[self.geomInd].max(axis=0)])
        else:
            extent = np.asarray(extent)  # Ensure np.array type.
            padGrid = False  # Skip padding option if grid preset.

        # Check for passed size - can be single or all dims.
        # Set to all dims if required.
        if extent.size < 3:
            extent = np.c_[np.ones(3)*extent[0], np.ones(3)*extent[1]].T

        # If set, force grid to square - same min and max for all dims
        if sqGrid:
            extent[0,:] = extent.min()
            extent[1,:] = extent.max()

        # Add grid padding
        if padGrid:
            extent[0,:] -= padGrid
            extent[1,:] += padGrid


        #*** Set axes
        # TODO: better storage options here - Xarray or HV?
        # See also efield.py for more.
        axes = []
        for col in extent.T:
            # print(col)
            axes.append(np.linspace(col[0], col[1], resolution))

        self.axes = np.array(axes)  # Set as axes here to avoid confusion with pyVista coords later.
        # self.x = axes[0]
        # self.y = axes[1]
        # self.z = axes[2]

        #*** Set grid
        # self.X, self.Y, self.Z = np.meshgrid(axes[0], axes[1], axes[2])  # Store to self
        X, Y, Z = np.meshgrid(axes[0], axes[1], axes[2])   # Local values only, used to set pyVista object
        # spacing = np.array((area_max - area_min)/resolution)

        # Set pyVista object to hold orbital calculations.
        # Note ordering!
        self.vol = pv.StructuredGrid(X, Z, Y)



    def calcOrb(self, orbN = None):
        """
        Calculate orbital volume.

        Parameters
        ----------
        orbN : int, optional, default = None
            Set orbital number to plot.
            This will default to the homo (as defined by cclib data.homos[-1]) if not set.
            Note there is a -1 indexing offset for zero-indexed python arrays, which is **applied** in the function.

        """

        if orbN is None:
            orbN = self.data.homos[-1]
        else:
            orbN -= 1 # Set to python indexing

        # Set orbital coeffs.
        # For unrestricted spin case, try stacking, athough this may be incorrect!
        if len(self.data.mocoeffs) == 1:
            coeffs = self.data.mocoeffs[0][orbN]
        else:
            print('***Warning: Urestricted orbitals not yet tested!')
            coeffs = np.array([self.data.mocoeffs[0][orbN], self.data.mocoeffs[1][orbN]])

        # Set function & calculate volume
        self.f[orbN+1] = self.wavefunction.molecular_orbital(self.data.atomcoords[self.geomInd], coeffs, self.data.gbasis)  # Note first input is atom positions, also in data.atomcoords

        # Calculate and assign as pyVista array, named as orbital N (note - must be str)
        # self.vol[orbN+1] = self.f(self.x, self.y, self.z).flatten(order='K')
        self.vol[str(orbN+1)] = self.f[orbN+1](self.vol.x, self.vol.y, self.vol.z).flatten(order='K')

        # Store in dictionary (or Xarray, or other...?)
        # AH - should just set pvObject here, and add each orbital as a data array.
        # Now done at init, and set by orbN above.


    def plotOrb(self, orbN = None, isoLevels = 6, isoValsAbs = None, isoValsPC = None, interactive = True, opacity = 0.5):
        """
        Plot orbital(s) with pyVista/ITK

        isoLevels : int, optional, default = 6
            Number of isosurfs to render.

        isoValsAbs : list or array, default = None
            Isovals to use (absolute values).

        isoValsPC : list or array, default = None
            Isovals to use, percentage values.

        """
        # For np.array types
        isoValsAbs = np.array(isoValsAbs)
        isoValsPC = np.array(isoValsPC)

        # Set ITK widgets or pv.Plotter
        if interactive:
            pl = pv.PlotterITK()
        else:
            pl = pv.Plotter()

        # Add atoms
        for n, item in enumerate(self.data.atomcoords[self.geomInd]):
        #     pl.add_mesh(pv.Sphere(center=item, radius=data.atomnos[n]/100))  # Relative scaling, small atoms
        #     pl.add_mesh(pv.Sphere(center=item, radius=np.sqrt(data.atomnos[n]/10)))  # sqrt scaling, space filling atoms
            pl.add_mesh(pv.Sphere(center=item, radius=np.sqrt(self.data.atomnos[n])/10))  # sqrt scaling, mid-sized atoms


        # Add meshes per orbital
        for item in self.vol.array_names:
            limitVals = [self.vol[item].min(), self.vol[item].max(), np.abs(self.vol[item]).mean()]

            if isoValsAbs is not None:
                isoValsOrb = np.array(isoValsAbs)
            elif isoValsPC:
                isoValsOrb = np.r_[isoValsPC, -isoValsPC] * limitVals[2]  # Set PC vals from mean?
            else:
                isoValsOrb = np.linspace(-limitVals[2], limitVals[2], isoLevels)  # Again set from mean.

            # Add contours for currently selected scalars (orbital)
            pl.add_mesh(self.vol.contour(isoValsOrb, scalars = item), smooth_shading=True, opacity=opacity)  # Plot iso = 0.1

        # Render plot
        self.pl = pl
        self.pl.show(True)
