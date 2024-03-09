"""
ePSproc electric field functions - polarization.

01/03/24 v1
    - Basic methods in place.
    - May be some ambiguities/issues in conversion between basis types.
        TODO: check phases more carefully.
    - Some functions require py_pol, but left as optional.
    - ep basis not yet tested, should be same as Elr except for linear case?
    - Only tested for single inputs so far, but should support multiples.


Quick pol overviews:

- https://en.wikipedia.org/wiki/Stokes_parameters#Representations_in_fixed_bases
- http://dielslab.unm.edu/sites/default/files/polarization_overview.pdf

py_pol library details: https://py-pol.readthedocs.io

TODO: implement additonal py_pol creation routines with wrapper?

"""

import numpy as np
import xarray as xr   # TODO: change to optional, since only used in setep

# Optional imports
py_polFlag = False
try:
    from py_pol.stokes import Stokes, create_Stokes, degrees
    py_polFlag = True

except ImportError as e:
    if e.msg != "No module named 'py_pol'":
        raise
    print('* py_pol not found, some field functions not available. \
          Run "pip install py_pol" to install.')


from epsproc.geomFunc import geomCalc
from epsproc.sphCalc import setPolGeoms, setBLMs, TKQarrayRotX
from epsproc.sphPlot import sphFromBLMPlot


class EfieldPol():
    """
    Class for handling E-field with polarization.

    Defined by either:

    - (Ex,Ey)
    - (El,Er)
    - ep list of field strengths as per geomFunc.EPR() routine.
        (NOT YET IMPLEMENTED for field creation, but use self.setep() for conversion.)

    If py_pol is installed, additional routines are available, including Stokes vector handling and plotting functions.

    Parameters
    ----------

    Exy : list or np.array, optional, default = None
        Field amplitudes defined by [Ex,Ey] pairs (list or array).
        Optionally can be in mag, phase format, with items [Ex,phix,Ey,phiy].
        And field will be set as [Ex*e^{i*phix}, Ey*e^{i*phiy}].

    Elr : list or np.array, optional, default = None
        Field amplitudes defined by [El,Er] pairs.
        E.g. Elr = [[El1,Er1],[El2,Er2]...]

    ell : list or np.array, optional, default = None
        Elliptical field defined as terms [[amplitude,azimuth,ellipticity]...]
        Only available if py_pol available, and then run Stokes.elliptical() method.
        NOTE azimuth in RAD.
        Ellipticity range: [-pi/4:pi/4]

    normField : bool, optional, default = True
        Renormalise input field components to unity magnitude if True.

    verbose : bool or int, optional, default = True
        Controls printing level to terminal.


    """

    # Init class with inputs for different methods.
    # May want to change to **kwargs to avoid overloaded locals() usage below.
    def __init__(self, Exy = None, Elr = None,
#                  ep = None, p = None,
                 ell = None,
                 normField = True, verbose = 1):

        # Set defaults
        self.inputs = locals()
        self.normField = normField
        self.verbose = verbose

        # Check for main config
        inputDict = {k:v for k,v in self.inputs.items() if v is not None and k not in ['self','normField','verbose']}

        if not inputDict:
            print(f"*** No inputs specified, setting default case Exy=[1,0]")
            inputDict['Exy'] = [1,0]

        # TODO: routine for setting class objs from dict here - see PEMtk.

        # Set methods per input
        # Assumes methods setEtype, and take one arg.
        # For more control set dict of definitions.
        for k,v in inputDict.items():
            func = getattr(self,f"set{k}",None)

            if func is not None:
                func(v)
                print(f"Set field from {k}.")

                if self.verbose:
                    print(getattr(self,k))
            else:
                print(f"*** Epol error: Function {k} not defined.")

#
#****** Basic functionality - Exy and Elr fields, py_pol wrappers
#

    def setExy(self, Exy):
        """Set Exy and dependent terms."""

        Exy = np.array(Exy, ndmin=2)  # Set ndmin here to allow singleton case
                                           # Assumes input Exy = [[Ex1,Ey1],[Ex2,Ey2]...]


        # For mag-phase definition
        if Exy.shape[1] == 4:
            self.Exy = np.c_[Exy[:,0]*np.exp(1j*Exy[:,1]),
                                 Exy[:,2]*np.exp(1j*Exy[:,3])]

        else:
            self.Exy = Exy

        # Norm field?
        # TODO: didn't test for complex case yet
        if self.normField:
#             self.norm = 1/np.sqrt((np.abs(self.Exy)**2).sum(axis=1))
#             self.Exy = self.Exy * np.tile(self.norm,(2,1)).T
            self.Exy, self.norm = self.calcNorm(self.Exy)

        # Set Elr
        # Defined as El = Ex - iEy, Er = Ex + iEy
#         norm = 1/np.sqrt(2)   # Always need extra norm here!
#         self.Elr = norm * np.c_[self.Exy[:,0]-1j*self.Exy[:,1],
#                                 self.Exy[:,0]+1j*self.Exy[:,1]]
        self.calcElrFromExy()

        if py_polFlag:
            self.stokes = Stokes('Light source from Exy')
            self.stokes.elliptical_light(self.Exy[:,0],self.Exy[:,1])


    def setElr(self,Elr):
        """
        Set Elr and dependent terms.

        Assumes input Elr = [[El1,Er1],[El2,E22]...]

        """

        Elr = np.array(Elr, ndmin=2)  # Set ndmin here to allow singleton case
                                      # Assumes input Elr = [[El1,Er1],[El2,Er2]...]

        # For mag-phase definition
        if Elr.shape[1] == 4:
            Elr = np.c_[Elr[:,0]*np.exp(1j*Elr[:,1]),
                        Elr[:,2]*np.exp(1j*Elr[:,3])]

        self.Elr = Elr

        # TODO: convert to Exy too
        # This is currently done for py_pol case only below.
        self.Exy = None

        if py_polFlag:
            # Set l+r components - may be a better way to do this?
            Sl = Stokes('LCP source')
            Sl.circular_light('l', intensity = self.Elr[:,0])
            Sr = Stokes('RCP source')
            Sr.circular_light('r', intensity = self.Elr[:,1])

            # Add components. Note global_phase, ensures El+Er for untiy case returns y-pol result.
            self.stokes = Sr+Sl.set_global_phase(180*degrees)

            # Also, if either component is empty, sum > NaNs unless phase is 0
            if np.isnan(self.stokes.M).any():
                self.stokes = Sr+Sl.set_global_phase(0*degrees)

#             if Sl.M.any() and Sr.M.any():
#                 self.stokes = Sr+Sl.set_global_phase(180*degrees)
#             elif Sl.M.any():
#                 self.stokes = Sl
#             else:
#                 self.stokes = Sr

            self.setExyFromStokes()


    def calcElrFromExy(self):
        # Set Elr
        # Defined as El = Ex - iEy, Er = Ex + iEy
        norm = 1/np.sqrt(2)   # Always need extra norm here!

        # Note sign convention here,
        # matches wiki defns. https://en.wikipedia.org/wiki/Stokes_parameters#Definitions_2
        self.Elr = norm * np.c_[self.Exy[:,0]+1j*self.Exy[:,1],
                                self.Exy[:,0]-1j*self.Exy[:,1]]


    def calcNorm(self,E):
        """
        Normalise E field

        Assumes 2D NP.array, e.g. Exy = [[Ex1,Ey1],[Ex2,Ey2]...], norm by row.

        """

        # Norm field?
        # TODO: didn't test for complex case yet
#         if self.normField:

        norm = 1/np.sqrt((np.abs(E)**2).sum(axis=1))
        E = E * np.tile(norm,(2,1)).T

        return E, norm


    def setell(self, elliptical):
        """
        Set elliptical light using py_pol

        elliptical = [[amplitude,azimuth,ellipticity]...]

        NOTE: currently assumes azimuth in RAD.

        From py_pol docs:

            azimuth (float): [0, pi]: azimuth. Default: 0
            ellipticity (float): [-pi/4, pi/4]: ellipticity. Default: 0
            intensity (numpy.array or float): Array of intensity. Default: 1.
            amplitude (numpy.array or float): Array of electric field amplitude. Overrides intensity if it is different than None. Default: None.

        See also https://py-pol.readthedocs.io/en/master/source/tutorial/Stokes.html#Parameters-of-Stokes-vector

        """

        if not py_polFlag:
            print("*** Please install py_pol to use `setelliptical()`. Run 'pip install py_pol'.")
            return None

        self.ell = np.array(elliptical,ndmin=2)

        S = Stokes('Elliptical light source')
        S.general_azimuth_ellipticity(amplitude=self.ell[:,0],
                                       azimuth=self.ell[:,1],   # *degrees,
                                       ellipticity=self.ell[:,2])

        self.stokes = S
        self.setExyFromStokes()
        self.calcElrFromExy()


    def setExyFromStokes(self):
        """Convert py_pol stokes params > Exy components."""
        # Set other basis terms
        # NOTE - force type to complex here, otherwise may accidentally be cast to real and drop phase term below!
        self.Exy = np.array(self.stokes.parameters.amplitudes(),ndmin=2).astype(np.complex)

        # Note transpose to give [[x,y]...] 2D array for multi field cases.
        if self.Exy.shape[0]>1:
            self.Exy = self.Exy.T

        # Propagate global phase term
        self.Exy[:,1] = self.Exy[:,1]*np.exp(1j*self.stokes.parameters.delay())


    def plot(self, **kwargs):
        if py_polFlag and hasattr(self,'stokes'):
            self.stokes.draw_ellipse(**kwargs)
        else:
            print("py_pol library required for self.plot. Run 'pip install py_pol' to install.")

#
#********** Additional functionality for ePSproc integration
#

    def setep(self, labels = None):
        """
        Set field terms in ep array format for use in :py:func:`epsproc.geomFunc.EPR()`

        Note this uses self.Elr ONLY, and sets p=[-1,1].
        Note EPR() currently only supports single pol state.

        self.epDict can be passed to EPR inputs.
        self.epXR contains Xarray representation.

        TODO: update EPR() for multiple pol states with labels.

        """

        p = [-1,1]

        self.epDict = {'ep':self.Elr.squeeze(),  # Note squeeze here for single pol state.
                       'p':p}


        if labels is None:
            # Set labels if missing
            labels = np.arange(1,self.Elr.shape[0]+1)

        self.epXR = xr.DataArray(self.Elr, coords={'p':p,'Labels':labels},
                                 dims=['Labels','p'])
#         self.epXR = xr.DataArray(self.Elr, coords={'p':p}, dims='p')

        if self.verbose:
            print("Set parameters to `self.epDict` and `self.epXR`.")


    def calcEPR(self):
        """
        Compute polarization tensor EPR from self.epDict.
        .. math:: E_{PR}(\hat{e})=[e\otimes e^*]^P_R

        Thin wrapper for :py:func:`epsproc.geomFunc.EPR()`

        """
        # Set inputs if missing.
        if not hasattr(self,'epDict'):
            self.setep()

        # Run tensor calculation.
        self.EPRX = geomCalc.EPR(form = 'xarray', **self.epDict)

        if self.verbose:
            print("Set geomCalc.EPR() results to `self.EPRX`.")


    def setOrientation(self, RX = None, eulerAngs = None, labels = None,
                        mapping = 'exy'):
        """
        Rotate Epol fields.

        In default case (no args), map Epol (x,y) to standard ePSproc axis system.

        Default mappings define 3 rotations:
              1. Epol(x) > z
              2. Epol(x) > x
              3. Epol(y) > y

        Where (1) is the usual z-axis definition for linear pol. case.
        And all 3 match the usual cases in ep.setPolGeoms().

        Default mappings are set as:

            pRot = [np.pi/2, 0, 0, np.pi/2, 0]
            tRot = [np.pi/2, 0, 0, np.pi/4, np.pi/4]
            cRot = [np.pi/2, np.pi/2, 0, np.pi/2, np.pi/2]
            labels = ['x>z', 'x', 'y', '+45y', '+45x']

            eulerAngs = np.array([pRot, tRot, cRot]).T

        """

        if RX is None:
            RX = setPolGeoms(eulerAngs = eulerAngs, labels = labels, defaultMap = mapping)

        self.RX = RX

        # Set terms - currently use BLM array for this.
        # TODO: tidy this up, test for multiple term passing.
        # TODO: pol array type?
        # Basic
        # self.YLM = setBLMs(np.array([[1,0,0],[1,-1,self.Elr[0,0]],[1,1,self.Elr[0,1]]]))
        # self.YLM = setBLMs(np.array([np.ones(self.Elr[:,0].size),self.Elr[:,0],self.Elr[:,1]]), LMLabels = np.array([[1,0],[1,-1],[1,1]]),
        # self.YLM = setBLMs(np.array([self.Elr[:,0],self.Elr[:,1]]), LMLabels = np.array([[1,-1],[1,1]]),
        #                     name = 'Epol')
        self.setYLM()

        # Rotate
        self.YLMrot, _, _ = TKQarrayRotX(self.YLM, RX)

        if self.verbose:
            print("Set pol state data to self.YLM and self.YLMrot, and orientations to self.RX.")

        # TODO: Set terms to standard dict - need to set options for output key in setep(), calcEPR()?


    def setYLM(self):
        """Set Elr terms to YLM expansion Xarray via :py:func:`epsproc.sphCalc.setBLMs`."""

        self.YLM = setBLMs(np.array([self.Elr[:,0],self.Elr[:,1]]), LMLabels = np.array([[1,-1],[1,1]]),
                            name = 'Epol')


    def plotSph(self, dataType = None, **kwargs):
        """
        Thin wrapper for :py:func:`epsproc.sphFromBLMPlot`

        Set working defaults. Note facetDim only working for Mpl backend currently (08/03/24).
        See base class padPlot() for better routine!

        Parameters
        ----------
        dataType : string, optional, default = None
            If None, plot from self.YLM
            If 'rot', plot from self.YLMrot

        """

        # Default case, note .squeeze() currently required, and not tested with multi-pol cases.
        if dataType is None:
            if not hasattr(self,'YLM'):
                self.setYLM()

            Itp, fig = sphFromBLMPlot(self.YLM.squeeze(drop=True), plotFlag = True, **kwargs)

        # For rotated frame plotted currently only working for facetDim='Euler'
        if dataType == 'rot':
            if not hasattr(self,'YLMrot'):
                print("Missing `self.YLMrot` data, run `self.setOrientation()` to define.")
            else:
                Itp, fig = sphFromBLMPlot(self.YLMrot.squeeze(drop=True), plotFlag = True,
                                    backend='mpl', facetDim='Euler', **kwargs)
