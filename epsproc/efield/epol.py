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

"""

import numpy as np

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
    
    
    
class EfieldPol():
    """
    Class for handling E-field with polarization.
    
    Defined by either:
    
    - (Ex,Ey)
    - (El,Er)
    - ep list of field strengths as per geomFunc.EPR() routine.
    
    If py_pol is installed, additional routines are available, including Stokes vector handling and plotting functions.
    
    Parameters
    ----------
    
    Exy : list or np.array, optional, default = None
        Field amplitudes defined by [Ex,Ey] pairs (list or array).
        Optionally can be in mag, phase format, with items [Ex,phix,Ey,phiy].
        
    Elr : list or np.array, optional, default = None
        Field amplitudes defined by [El,Er] pairs.
        E.g. Elr = [[El1,Er1],[El2,Er2]...]
        
    elliptical : list or np.array, optional, default = None
        Field defined as terms [[amplitude,azimuth,ellipticity]...]
        Only available if py_pol available, and then run Stokes.elliptical() method.
        NOTE azimuth in degrees.
        
    
    """
    
    def __init__(self, Exy = None, Elr = None, 
#                  ep = None, p = None,
                 elliptical = None,
                 normField = True):
        
        # Set defaults
        self.inputs = locals()
        self.normField = normField
        
        # Check for main config
        inputDict = {k:v for k,v in self.inputs.items() if v is not None and k not in ['self','normField']}
        
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
                print(f"Set field {k}.")
            else:
                print(f"*** Epol error: Function {k} not defined.")
        
        
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
        
        self.Elr = np.array(Elr, ndmin=2)  # Set ndmin here to allow singleton case
                                      # Assumes input Elr = [[El1,Er1],[El2,Er2]...]
        
        # TODO: convert to Exy too
        # This is currently done for py_pol case only below.
        
        
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
        self.Elr = norm * np.c_[self.Exy[:,0]-1j*self.Exy[:,1],
                                self.Exy[:,0]+1j*self.Exy[:,1]] 
    
    
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
    
    
    def setelliptical(self, elliptical):
        """
        Set elliptical light using py_pol
        
        elliptical = [[amplitude,azimuth,ellipticity]...]
        
        NOTE: currently assumes azimuth in DEGREES.
        
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
                                       azimuth=self.ell[:,1]*degrees, 
                                       ellipticity=self.ell[:,2])
 
        self.stokes = S
        self.setExyFromStokes()
        self.calcElrFromExy()
        
    
    def setExyFromStokes(self):
        """Convert py_pol stokes params > Exy components."""
        # Set other basis terms
        self.Exy = np.array(self.stokes.parameters.amplitudes(),ndmin=2)
        
        # Propagate global phase term
        self.Exy[:,1] = self.Exy[:,1]*np.exp(1j*self.stokes.parameters.delay())
        
        
    def plot(self):
        if py_polFlag and hasattr(self,'stokes'):
            self.stokes.draw_ellipse()
        else:
            print("py_pol library required for self.plot. Run 'pip install py_pol' to install.")