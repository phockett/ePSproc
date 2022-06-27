"""
ePSproc Efield class

"""

import numpy as np
from scipy import constants as scipy_constants
import pprint

import xarray as xr

# Additional imports for plotting
import matplotlib.pyplot as plt

import holoviews as hv
from holoviews import opts



# Use dictionary to create field defns intially?  Gets confusing otherwise.
# Should be able to iterate/unpack values here for different cases.
# Minimal defn. with (sigma, wavelength) or (bandwidth, wavelength), in working units.
# Edef = {'sigma': 15,
#         'l0': 800
#         }


class Efield():
    '''
    Basic class for handling E-field generation, data and related functions.

    Currently set for field creation based on Gaussian defined in time-domain. Defaults to a single (E,t) point if no parameters are supplied - this can be used with :py:func:`epsproc.geomFunc.EPR`

    Starting to look a bit top-heavy... passing directly to Xarray might save some structural overhead here?

    Define Gaussian field in the time-domain, with amplitude $E_0$, width $\sigma$ and carrier frequency $\omega_0$:

.. math::
    E(\omega_0,t)=E_{0}e^{-(t^2/2\sigma^{2})}e^{i\omega_0 t}

    Carrier frequency in rad/s:

.. math::
    \omega_0 = 2\pi \Omega_0

    FWHM and Gaussian width (see https://en.wikipedia.org/wiki/Gaussian_function for full definitions):

.. math::
    FWHM = 2\sigma\sqrt{2ln2}


    Parameters
    ----------
    If no parameters are passed, return single point with E0 = 1
    For a time-dependent field, pass at least sigma.

    E0 : float, optional, default = None
        Set (relative) field strength.
        If None, set E0 = 1

    sigma : float, optional, default = None
        Set Gaussian width :math:\sigma, where :math:FWHM=2\sigma\sqrt{2ln2}

    t : array, optional, default = None
        Set time axis for E-field.
        If None, defaults to np.arange(-(5*sigma), 5*sigma, self.dt), and self.dt = 0.01*sigma

    dt : float, optional, default = None
        Set time-step. Default will set from t-axis, or self.dt = 0.01*sigma if sigma only passed.
        Pass here to override.

    l0 : float, optional, default = None
        Central wavelength in wavelength working units.
        If supplied this will be used to set f0 and w0.

    f0, w0 : float, optional, default = None
        Carrier frequency in freq. units, and :math:\omega_0 in (rad/s).
        If only f0 supplied, w0 will be calculated. If f0 and w0 supplied, w0 will be used directly.
        If None, this will be ignored.

    units : dictionary, optional, default = None
        Set units. If not passed, default dict will be created with units (THz, nm, fs) set.


    TODO:
    - best way to handle unit conversions?  Did some of this for FROG code already? Yes, see python\test_codes\E_fields_redux_231118_class.py
    - Consider other structures here, dataclass or namedtuple? https://stackoverflow.com/questions/354883/how-do-i-return-multiple-values-from-a-function
    - Xarray conversion?  Use Xarray directly in class?
    - Bug with FFT/iFFT code?  If shift=True and ishift=True then get extra freq components on iFFT - must be bug in axis assignment/conversion?  OK if shift = False.

    '''

    def __init__(self, Edef = None, units = None, verbose = True):  # , E0 = None, sigma = None, A = 1, t = None, dt = None, l0 = None, f0 = None, w0 = None, units = None):

        #*** Define master list of parameters & init as empty or with defaults
        # If no parameters are passed, this can be used to create a defn template
        self.Edef = {'Pulse' : {'type' : 'Gaussian',   # Eventually should generalise this to allow for other pulse types.
                                'domain' : 't',       # Set a domain label here to switch between (t,E) domain pulse defn. - MAY WANT TO MOVE THIS TO SEPARATE LIST FN?
                                'dfft' : 'f',        # Set conjugate domain (for FFT)
                                'sigma' : None,
                                'E0' : 1,       # Set magnitude of field
                                'CEP':0,        # Set CEP, only used if pulse carrier freq is also defined
                                'A' : 1,        # Set A for Gaussian - should just set as E0?
                                'FWHM' : None,
                                'p' : 0    # Set polarization state. For single values {-1,0,+1} this will just label output. For 'XY' calc. spherical terms.
#                                 'origin' : 0    # Set origin for domain
                               },

                     'Freq' : {'l' : None,   # Set defns. for carrier, (wavelength, freq, ang. freq). Assumed to be in working units. Also use to index units later.
                               'f' : None,
                               'w' : None,
                               # 'CEP':0,        # Set CEP, only used if pulse carrier freq is also defined - SET IN PULSE for now, due to use of assignment loop for Freq terms.
                               # 'Ehv' : None  # Should add here, but need to fix assignment routine first!  Now set by fConv() method.
                               },

                     'Ehv' : {'Ef' : None,  # Defined energy domain to hold field converted to E(hv) units
                              'axis' : None
                             },

                     't' : {'Ef' : None,
                            'axis' : None,    # Set generic labels here to allow for (t,f) axis shared methods.
                            'delta' : None
                               },

                     'f' : {'Ef' : None,
                            'axis' : None,
                            'delta' : None
                               },

#                      'EField' : {'Et' : None,
#                                  'Ef' : None
#                                  },

#                      'Spectrogram' : {'gate' : None,
#                                       'data' : None
#                                      },

                     'Units' : {   # 'f':{'value':1e12, 'unit':'THz'},  # Set this as a derived unit, 1/t, below, otherwise may get discrepancies in unit conversion.
                                'l':{'value':1e-9, 'unit':'nm'},
                                't':{'value':1e-15, 'unit':'fs'},
                                'Ehv':{'value':1, 'unit':'eV'}
                                },

                     # Global FFT settings
                     # May have issues here with shift and ishift - always need the latter...?
                     # Bug with FFT/iFFT code?  If shift=True and ishift=True then get extra freq components on iFFT - must be bug in axis assignment/conversion?  OK if shift = False.
                     'FFT' : {'shift':False, 'ishift':True,
                              'pad':True, 'positiveHalf':False,
                              'phaseMaskFlag':False, 'thres':1e-3}
                    }

        self.Emod = {'N':0, 0:{}}  # Set empty dict to hold modified fields (propagated, shaped etc.). Use N to hold next field to update value.


        # Assign any passed values
        # Skip if no Edef is passed (will just set blank dict)
        if Edef is not None:
            for key in Edef:
                if key in self.Edef:
                    for subkey in Edef[key]:
                        if subkey in self.Edef[key]:
#                             print(f'Setting Item {subkey} in dict {key} to {Edef[key][subkey]}')
                            self.Edef[key][subkey] = Edef[key][subkey]
                        else:
                            print(f'Item {subkey} in dict {key} not recognised.')
                else:
                    print(f'Key {key} not recognised.')

            #*** Set working units for conversions
            # Set f units
            self.Edef['Units']['f'] = {'value': np.round(1/self.Edef['Units']['t']['value']), 'unit': f"1/{self.Edef['Units']['t']['unit']}"}  # Derived unit

            # Set as Hz in specific cases
            if self.Edef['Units']['f']['value'] == 1e15:
                self.Edef['Units']['f']['unit'] = 'PHz'
            if self.Edef['Units']['f']['value'] == 1e12:
                self.Edef['Units']['f']['unit'] = 'THz'

            # Set units for E - NOW SET DERIVED UNITS via fConv() method
#             self.Edef['Units']['E'] = self.Edef['Units']['f'].copy()

            # Set c in working units
            self.Edef['Units']['c'] = {'value': scipy_constants.c * self.Edef['Units']['t']['value']/self.Edef['Units']['l']['value'],
                                       'unit': f"{self.Edef['Units']['l']['unit']}/{self.Edef['Units']['t']['unit']}"}
#             self.Edef['Units']['c']['value'] = scipy_constants.c * self.Edef['Units']['t']['value'] /self.Edef['Units']['l']['value']

            #**** Calculate field and derived properties if set
            self.setEf()

            # Set description string - use """ format string, or ""\ multiline
            domain = self.Edef['Pulse']['domain']
            self.Estring = f"""{self.Edef['Pulse']['type']} pulse: $\sigma$={self.Edef['Pulse']['sigma']} {self.Edef['Units'][domain]['unit']}, FWHM={self.Edef['Pulse']['FWHM']:.3f} {self.Edef['Units'][domain]['unit']},
l0={self.Edef['Freq']['l']:.3f} (dl={self.Edef['Pulse']['dl']:.3f}) {self.Edef['Units']['l']['unit']},
f0={self.Edef['Freq']['f']:.3f} (df={self.Edef['Pulse']['dw']:.3f}) {self.Edef['Units']['f']['unit']}
(bandwidths for Gaussian transform limited pulses)"""

            # Print summary details
            if verbose:
                self.printDef()


    def printDef(self):
        print('Pulse properties dictionary set:')
        pp = pprint.PrettyPrinter(indent=4)

        # Print defns.
        for key in ['Pulse', 'Freq', 'Units', 'FFT']:
            print(f'{key} settings:')
            pp.pprint(self.Edef[key])

        # Print domain details
        for key in [self.Edef['Pulse']['domain'], self.Edef['Pulse']['dfft']]:
            if self.Edef[key]['Ef'] is not None:
                print(f'{key} domain settings:')
                print(f"Points = {len(self.Edef[key]['Ef'])}, delta = {self.Edef[key]['delta']}")
            else:
                print(f"{key} domain not set")

        # Print field details
        for key in ['EField', 'Spectrogram']:
            pass


    def setEf(self):
        """
        Calculate pulse properties and E-fields based on Edef.
        """

        # Set pulse in time-domain
        if self.Edef['Pulse']['domain'] == 't':

            # Set pulse FWHM or sigma if defined.
            # If both values are preset, only sigma is used.
            if (self.Edef['Pulse']['sigma'] is not None) and (self.Edef['Pulse']['type'] == 'Gaussian'):
                self.FWHMGau()
            elif(self.Edef['Pulse']['FWHM'] is not None) and (self.Edef['Pulse']['type'] == 'Gaussian'):
                self.sigmaGau()

            # Set t axis if required
            if (self.Edef['t']['axis'] is None) and (self.Edef['Pulse']['sigma'] is not None):

                if self.Edef['t']['delta'] is None:
                    self.Edef['t']['delta'] = 1e-3*self.Edef['Pulse']['sigma']  # Default step size for t axis, relative to sigma. May get large however!

                self.Edef['t']['axis'] = np.arange(-(5*self.Edef['Pulse']['sigma']), 5*self.Edef['Pulse']['sigma'], self.Edef['t']['delta'])

            # No pulse in this case, just set to single point (t=0)
            elif (self.Edef['t']['axis'] is None) and (self.Edef['Pulse']['sigma'] is None):
                self.Edef['t']['axis'] = [0]

            # If t-axis is passed, set dt from this - assumes linear axis
            # Use length here to allow for list or np.array types
            if (len(self.Edef['t']['axis']) > 1) and (self.Edef['t']['delta'] is None):
                self.Edef['t']['delta'] = np.abs(self.Edef['t']['axis'][1]-self.Edef['t']['axis'][0])


            # Check and set carrier freq if set
            if any(self.Edef['Freq'].values()):
                # Check which value is set, and set missing values.
                # Neater/slicker way to do this? Recursively?
                for key in self.Edef['Freq']:
                    if self.Edef['Freq'][key] is not None:
                        refKey = key
                        refValue = self.Edef['Freq'][key]
                        # print(refValue)

                if refKey == 'w':
                    refValue *= self.Edef['Units']['f']['value']/2*np.pi

# TODO: set for Ehv case.
#                 if refKey is 'Ehv':
#                     refValue *= self.Edef['Units']['f']['value']/2*np.pi

#                 else:

#                 refValue *= self.Edef['Units'][refKey]['value']  # Convert to real units - not required if c in working units

                for key in self.Edef['Freq']:
                    if key not in [refKey, 'w']:
#                         self.Edef['Freq'][key] = (scipy_constants.c/refValue)/self.Edef['Units'][key]['value']  # Convert between wavelength and freq., and set to working units
                        self.Edef['Freq'][key] = self.Edef['Units']['c']['value']/refValue  # With c in working units

                if refKey != 'w':  # Set w0 in working units (rad/[unit t])
#                     self.Edef['Freq']['w'] = 2*np.pi * self.Edef['Freq']['f'] * self.Edef['Units']['f']['value'] * self.Edef['Units']['t']['value']
                    self.Edef['Freq']['w'] = 2*np.pi * self.Edef['Freq']['f']


            self.ECalc()   # Set defined field
            self.EFFT()    # Set FFT field

            # Set fields in Emod too - may eventually replace above with this?
            for key in [self.Edef['Pulse']['domain'], self.Edef['Pulse']['dfft']]:
                self.Emod[self.Emod['N']][key] = self.Edef[key]

            self.Emod['N'] += 1



    #***************** Basic generators

    # Define Gaussian pulse in time domain, if carrier freq. is not defined calc. envelope only
    # Sigma = Gaussian width, FWHM = 2
    # Define Gaussian
    def Gau(self):
#             g = np.exp(-0.5*self.Edef['Pulse']['A']*(self.Edef['t']['axis']**2)/(self.Edef['Pulse']['sigma']**2))
#         g = np.exp(-0.5*self.Edef['Pulse']['A']*((self.Edef[self.Edef['Pulse']['domain']]['axis']/self.Edef['Pulse']['sigma'])**2))
        g = np.exp(-0.5*self.Edef['Pulse']['A']*(self.Edef[self.Edef['Pulse']['domain']]['axis']**2)/(self.Edef['Pulse']['sigma']**2))

        return g

    # Define FWHM from sigma, for a Gaussian pulse
    def FWHMGau(self):
        self.Edef['Pulse']['FWHM'] = 2*np.sqrt(2*np.log(2))*self.Edef['Pulse']['sigma']
        self.Edef['Pulse']['dw'] = 0.44/self.Edef['Pulse']['sigma']  # Spectral width for a Gaussian pulse, working units, sigma*tUnit to give in Hz
                                                                     # This is a bit of a fudge - need to decide where to put dw, and E domain defn.

        if self.Edef['Freq']['l'] is not None:
            self.Edef['Pulse']['dl'] = (self.Edef['Pulse']['dw']*self.Edef['Freq']['l']**2)/self.Edef['Units']['c']['value']
#             self.Edef['Pulse']['dl'] = ((self.Edef['Pulse']['dw']*self.Edef['Freq']['l']**2)/scipy_constants.c)  # Spectral width in wavelength units
#         print(f"For sigma={self.Edef['Pulse']['sigma']}: FWHM={self.Edef['Pulse']['FWHM']:.3f}, spectral width (transform limit)={self.Edef['Pulse']['dw']:.3f}")

    # Set field based on defined pulse parameters
    def ECalc(self):

        domain = self.Edef['Pulse']['domain']

        self.Edef[domain]['Ef'] = self.Edef['Pulse']['E0'] * self.Gau()

        if (self.Edef['Freq']['w'] is not None) and (domain == 't'):  # Only valid for t-domain defn.
#             print('OK')
#             print(self.Edef['Freq']['w'])
#             print(self.Edef[domain]['axis'])
#             print(np.exp(1.0j*self.Edef['Freq']['w']*self.Edef[domain]['axis']))
#             self.Edef[domain]['Ef'] = self.Edef[domain]['Ef'] * np.exp(1.0j*self.Edef['Freq']['w']*self.Edef[domain]['axis']/self.Edef['Units'][domain]['value'])
#             self.Edef[domain]['Ef'] = self.Edef[domain]['Ef'] * np.exp(1.0j*self.Edef['Freq']['w']*self.Edef[domain]['axis']*self.Edef['Units'][domain]['value'])
            self.Edef[domain]['Ef'] = self.Edef[domain]['Ef'] * np.exp(1.0j*(self.Edef['Freq']['w']*self.Edef[domain]['axis'] - self.Edef['Pulse']['CEP']))


    #********************* FT functions

    # Set other domain field as FFT
    # Define spectral domain pulse as FFT(E(t)). Return normalised freq. axis if t and dt are defined.
    def EFFT(self): #, pad = True, positiveHalf = True, phaseMaskFlag = False, thres = 1e-3):  # Now set in dict.

        domain = self.Edef['Pulse']['domain']
        dfft = self.Edef['Pulse']['dfft']

        # Use zero-padding for higher resolution FFT result?
        if self.Edef['FFT']['pad'] is True:
            n = np.int(2**(np.ceil(np.log2(self.Edef[domain]['Ef'].shape[0]))+3))
            nAxis = n
        else:
            n = None # Default value to pass to np.fft for no padding
            nAxis = self.Edef[domain]['axis'].shape[-1]

        Ebar = np.fft.fft(self.Edef[domain]['Ef'], n=n)   # No fft shift, axis [0 .... +ve .... -ve], see https://numpy.org/doc/stable/reference/generated/numpy.fft.fftfreq.html?highlight=fft
        axis = np.fft.fftfreq(nAxis, d=self.Edef[domain]['delta'])

        if self.Edef['FFT']['shift']:
            Ebar = np.fft.fftshift(Ebar)   # Apply fft shift to move 0 to centre of range
            axis = np.fft.fftshift(axis)

        # Set for full or half FT
        if self.Edef['FFT']['positiveHalf']:
            if self.Edef['FFT']['shift']:
                inds = np.arange(np.int(Ebar.shape[0]/2), Ebar.shape[0])   # Set for half the range, starting from centre
            else:
                inds = np.arange(0, np.int(Ebar.shape[0]/2))   # Set for half the range, starting at 0

            self.Edef[dfft]['Ef'] = Ebar[inds]
            self.Edef[dfft]['axis'] = axis[inds]
            self.Edef[dfft]['delta'] = axis[1]-axis[0]

        else:
            self.Edef[dfft]['Ef'] = Ebar
            self.Edef[dfft]['axis'] = axis
            self.Edef[dfft]['delta'] = axis[1]-axis[0]

        if self.Edef['FFT']['phaseMaskFlag']:
            self.phaseMask(domain = dfft)


#     # Function for checking for modified fields & sorting
#     def checkEmod(self):
# #         # Loop over domains
# #         for domain in [self.Edef['Pulse']['domain'], self.Edef['Pulse']['dfft']]:

# #             # Check for existing mod fields
# #             if domain in self.Emod.keys():
# #                 N = list(self.Emod[domain].keys())

# #         # Should change formatting here, may want domain_N to allow for multiple cases, and looping.
# #         if (domain + '_mod') in self.Edef.keys():
# #             N = list(self.Edef[domain + '_mod'].keys())
# #         else:
# #             self.Edef[domain + '_mod'] = {}
# #             N = 0

#         # Check if entries exist
#         if self.Emod:
#             N = self.Emod.N
#         else:
#             N = 0

#         # Set

    # Calculate iFFT(E).
    # This, sort of, assumes that the field is defined as E(w), and the spectral phase is modified.
    # But could be used for other cases.
    # Emod is used if passed, otherwise uses self.Edef[dfft]['Ef']
    # 27/04/20 Modified to use self.Emod for fields.
    # In this case, send either Emod to use this directly, or Ninput to use self.Emod[Ninput].
    # Set comment = '' to pass notes on field generation.
    def EiFFT(self, Emod = None, Ninput = None, comment = ''): #, f=None, pad=False):

        domain = self.Edef['Pulse']['domain']
        dfft = self.Edef['Pulse']['dfft']

        Nlast = self.Emod['N'] - 1  # Most recent field index, use this as input if nothing else specified

        # Set field based on input - this is currently a bit ugly!
#         if (Emod is None) and (Ninput is None):
#             if dfft in self.Emod[Nlast]:
#                 Emod = self.Emod[Nlast][dfft]['Ef']   # Default to most recent field if it exists
#                 EmodAxis = self.Emod[Nlast][dfft]['axis']
#             else:
#                 Emod = self.Edef[dfft]['Ef']  # Revert to original field defn. if not supplied
#                 EmodAxis = self.Edef[dfft]['axis']

#         elif (Emod is None) and (Ninput is not None):  # Set specific field from Emod
#             Emod = self.Emod[Ninput][dfft]['Ef']
#             EmodAxis = self.Emod[Ninput][dfft]['axis']

        # Rewrite for updated Emod dict.
        N = Nlast  # Set default
        if Ninput is not None:
            N = Ninput

        if Emod is None:
            Emod = self.Emod[N][dfft]['Ef']   # Default to most recent field if it exists
            EmodAxis = self.Emod[N][dfft]['axis']

            if 'comment' in self.Emod[N][dfft].keys():
                comment += self.Emod[N][dfft]['comment']  # Propagate comment

        # Transform back to time-domain, and center
#         Eifft = np.fft.ifftshift(np.fft.ifft(Emod))   # With shift
        Eifft = np.fft.ifft(Emod)  # Without shift
        axis = np.fft.fftfreq(Eifft.shape[-1], d=self.Edef[dfft]['delta'])

        if self.Edef['FFT']['ishift']:
            Eifft = np.fft.ifftshift(Eifft)  # Apply fft shift if set
            axis = np.fft.ifftshift(axis)

        # Set for full or half FT - if used with iFFT shift and a pulse center at 0, this will slice result.
        inds = np.arange(0, Eifft.shape[0])   # Set for full range

        if self.Edef['FFT']['positiveHalf']:
            Eifft /= 2.0   # Correct amplitudes, get x2 otherwise if positiveHalf is True.

            if self.Edef['FFT']['ishift']:
                inds = np.arange(np.int(Eifft.shape[0]/2), Eifft.shape[0])   # Set for half the range, starting from centre
            else:
                inds = np.arange(0, np.int(Eifft.shape[0]/2))   # Set for half the range, starting at 0
#         self.Edef[domain]['Ef'] = np.c_[self.Edef[domain]['Ef'], Eifft]   # Basic stacking OK if pad=False and positiveHalf=False, otherwise axis lengths different


#         if self.Edef['FFT']['positiveHalf']:
#             inds = np.int(Eifft.shape[0]/2)   # Set for half the range

#         # Set as new field, or as new dict?
#         # Should change formatting here, may want domain_N to allow for multiple cases, and looping.
#         if (domain + '_mod') in self.Edef.keys():
#             N = list(self.Edef[domain + '_mod'].keys())
#         else:
#             self.Edef[domain + '_mod'] = {}
#             N = 0

        # Set pair of fields in output
        # Always set to new output, or check for pair...?
        if domain in self.Emod[Nlast]:  # If conjugate domain is already set, create a new field pair...
            Noutput = Nlast + 1
            self.Emod[Noutput] = {}
            self.Emod[Noutput][dfft] = {'Ef':Emod,
                                        'axis':EmodAxis,
                                        'domain':domain
                                        }
        else:
            Noutput = Nlast  # ...otherwise use input index.

        self.Emod[Noutput][domain] = {'Ef':Eifft[inds],
                                      'axis':axis[inds],
                                      'domain':domain,
                                      'comment':comment  # Add passed comment here, may also want to autogenerate note on field source?

                                    }



        #if f:
        #    t = np.fft.ifftfreq(f.shape[-1])
        #else:
        #    t = None

    # TODO
    def finst(self):
        """Calculate f(t) as instantaneous freq."""
        print('Not impemented')


    #*************** Phase modification fns.
    # TODO: add more sophisticated fns, see, e.g., froglib.phasemanipulations, for removing linear and phase offsets.

    # Mask phase away from pulse?
    def phaseMask(self, domain = None, thres = None):
        """Mask pulse phase away from main features.  Should convert to use masked arrays for reversibility, or just return a mask (currently overwrites Ef)."""

        if thres is None:
            thres = self.Edef['FFT']['thres']  # Use global thres if not passed.

        thres = thres * np.abs(self.Edef[domain]['Ef']).max()

        self.Edef[domain]['Ef'][np.abs(self.Edef[domain]['Ef']) < thres] = 0  # Basic mask to zero away from features
#         self.Edef[domain]['Ef'][np.abs(self.Edef[domain]['Ef']) < thres] = np.nan  # Set as NaNs - will remove points from plots, but also throws errors with np.angle()


    def setPhase(self, phaseVec = None, domain = None):
        """Set imaginary field terms (set phase) for specified domain & recalculated FFT."""

        # Example: set quadratic phase
        # Set quadratic phase (corresponds to linear chirp), phi=A*f^2
        # A = 0.1/fUnit
        # Ewq = np.abs(E2w)*np.exp(1.0j*(A*(f-(f0/tUnit)))**2)

        # Transform back to time-domain
        # E2tmod = np.fft.ifftshift(np.fft.ifft(Ewq))

        pass
        # TODO: impement, and add original & modified field placeholders. Or set as new rows in Ef array?



    def removePhase(self, domain = None):
        """Remove imaginary field terms (reset phase)."""

        if domain is None:
            for item in ['t', 'f']:
                self.Edef[item]['Ef'] = np.abs(self.Edef[item]['Ef'])
        else:
            self.Edef[domain]['Ef'] = np.abs(self.Edef[domain]['Ef'])



    def chirp(self, A, resetPhase=False, comment = None):
        """
        Add quadratic phase (chirp) in spectral domain. Requires an E-field object, and chirp parameter A.

        .. :math: phi=A*(f-f0)^2

        TODO: check Diels for nomenclature here.

        """
        domain = 'f'

        Nlast = self.Emod['N']  # Last output slot, use N-1 as input if nothing else specified

        # Set field - checks no longer required as now Emod[0] set at init (in setEf()).
#         try:
#             Ew = self.Emod[Nlast-1][domain]['Ef']
#             EwAxis = self.Emod[Nlast-1][domain]['axis']
#         except KeyError:
#             Ew = self.Edef[domain]['Ef']
#             EwAxis = self.Edef[domain]['axis']

        Ew = self.Emod[Nlast-1][domain]['Ef']
        EwAxis = self.Emod[Nlast-1][domain]['axis']


        # Remove existing phase - this is only the linear term if starting from a Gaussian pulse
        if resetPhase:
            Ew = np.abs(Ew)

        # Add chirp
#         Ewq = np.abs(Ew)*np.exp(1.0j*(A*(self.Edef[domain]['axis']-self.Edef['Freq']['f']))**2)
#         Ewq = Ew*np.exp(1.0j*A*(self.Edef[domain]['axis']-self.Edef['Freq']['f'])**2)
        Ewq = Ew*np.exp(1.0j*A*(EwAxis-self.Edef['Freq']['f'])**2)


        # Set as new field, or as new dict?
        # Should change formatting here, may want domain_N to allow for multiple cases, and looping.
        # Or nest this? Or concatenate? (-- Only if domain axes are identical.)
#         if (domain + '_mod') in self.Edef.keys():
#             N = list(self.Edef[domain + '_mod'].keys())
#         else:
#             self.Edef[domain + '_mod'] = {}
#             N = 0

#         self.Edef[domain + '_mod'][N] = {'Ef':Ewq,
#                                      'axis':self.Edef[domain]['axis'],
#                                      'domain':domain

#                                     }

        self.Emod[Nlast] = {}
        self.Emod[Nlast][domain] = {'Ef':Ewq,
                                     'axis':EwAxis,
                                     'domain':domain,
                                     'comment': f'$E_{{{Nlast-1}}}$, phase {domain} chirped, A={A}',
                                     'A':A
                                    }

        self.Emod['N'] += 1  # Update indexer

#         self.EiFFT(Emod = Ewq)  # Set tmod pulse via ifft
        self.EiFFT()  # Set tmod pulse via ifft


    #***************** Spectrograms

    # Basic (Frog) spectrograms - code adapted from Froglib, https://github.com/xmhk/froglib
    # For code implementing various Frog methods, see "E_fields_redux_231118_class.py" - to be implemented here as frog() method
    def calcSpectrogram(self, signal = None, gate = None, fType = 'blind', N = None, verbose = False):

        domain = 't'

        # Set signal and gate fields. If not passed, use most recently set fields.
        if N is None:
            N = self.Emod['N'] - 1

        if signal is None:
            signal = self.Emod[N][domain]['Ef']

        elif type(signal) is np.ndarray:
            pass

#         elif type(signal) is   ####### May want to allow for passing of pulse defn. dictionary here?

        gateObj = None

        if gate is None:
            gate = signal   # Default to signal, or to short gaussian...?

        elif type(gate) is float:   # Take a passed value as a Gaussian width...?
            pass

        elif type(gate) is dict:  # Generate gate pulse as new Ef object
            gateObj = Efield(gate, verbose = verbose)
            gate = gateObj.Edef[domain]['Ef']

#         elif type(signal) is np.ndarray:
#             pass

        # TODO:
        # - Error checking, currently needs square array.
        # - Downsampling for cases with large FFT axis - just select ROI around features (see Frog code?)
        # - Methods, see "E_fields_redux_231118_class.py" for more frog types.

        # Following code in froglib...
        nn = len(signal)
        n2 = int(nn / 2)

        # (1) Outer product of two fields (time-domain), blind Frog case
#         ap = np.outer(signal, gate)

        # (1) Outer product of two fields (time-domain), depending on type of Frog
        # NOTE field ordering matters for X-Frog definitions.
        # Set options using dictonary (see https://simonwillison.net/2004/May/7/switch/ and https://stackoverflow.com/questions/60208/replacements-for-switch-statement-in-python)
        ap = {
              'SHG':    lambda f1,f2: np.outer(f1, f2) + np.outer(f2, f1),  # SHG case, symmetric in time
              'blind':  lambda f1,f2: np.outer(f1, f2),                     # Blind case, just two fields
              'SD':     lambda f1,f2: np.outer(f1**2,f2.conjugate()),       # SD classic
              'SDr':    lambda f1,f2: np.outer(f2**2,f1.conjugate()),       # SD classic - field ordering reversed, matters in X-Frog case
        #      'SDb':    lambda f1,f2: np.outer(f1,f2.conjugate()**2),
              'SD1':    lambda f1,f2: np.outer(f1, f1.conjugate()*f2),      # Various options for alternative SD cases, depending on non-linearity and field ordering
              'SD2':    lambda f1,f2: np.outer(f1*f2,f2),
              'SD3':    lambda f1,f2: np.outer(f1,f1*f2.conjugate()),
              'PG':     lambda f1,f2: np.outer(f1, np.abs(f2)**2),          # PG classic - results match Trebino for cubic case (flipped relative to SD)
              'PGr':    lambda f1,f2: np.outer(f2, np.abs(f1)**2),         # PG classic - field ordering reversed, matters in X-Frog case
              'PG1':    lambda f1,f2: np.outer(np.abs(f1)**2, f2),          # PG classic - this defn. identical to SD case.
              'TG1':    lambda f1,f2: np.outer(f1, f2**2)                   # TG options
              # 'TG2':    lambda f1,f2: np.outer(f1, f2**2)
              }[fType](signal, gate)

        # (2) Defined empty arrays to hold results
        m1 = np.zeros(np.shape(ap), dtype=np.complex128)
        m2 = np.zeros(np.shape(ap), dtype=np.complex128)

        # (3) Loop over input and roll - effectively sets tau for each row
        for i in range(n2 - 1, -n2, -1):
            m1[i + n2, :] = np.roll(ap[i + n2, :], -i)

        m1 = np.transpose(m1)

        # (4) Roll and FFT to set a freq. axis
        for i in range(nn):
            m2[i, :] = np.roll(np.fft.fft(np.roll(m1[i, :], +n2)), -n2)  # time-freq

        m2 = np.transpose(m2)  # freq - time
        m2 = m2 / np.max(np.max(np.abs(m2)))

#         return m2

        # Set outputs - should just set in Emod....?
        self.Spectrogram = {'siganl':signal,
                            'gate':gate,
                            'gateObj':gateObj,
                            'data':m2,
                            'N':N
                           }


    #***************** Derived domains/unit conversion

    # Set other domains via copy & rescale - would be neater just to store multiple axes, but duplicate for now.
    def fConv(self):
        """Convert freq. domain axis to lenght & energy units. and set fields."""

        domain = 'f'

        for N in np.arange(0, self.Emod['N']):

            # Set wavelength scale
            self.Emod[N]['l'] = self.Emod[N][domain].copy()   # Without .copy() this will just be a pointer.
            self.Emod[N]['l']['axis'] = self.Edef['Units']['c']['value']/self.Emod[N]['f']['axis']

            # Set energy scale
            self.Emod[N]['Ehv'] = self.Emod[N][domain].copy()
            self.Emod[N]['Ehv']['axis'] = (scipy_constants.h/scipy_constants.e) * self.Emod[N]['f']['axis']/self.Edef['Units']['t']['value']



    #***************** Conversion...

    # Convert a single set of fields to an Xarray Dataset
    # With spectrogram + looping over domains - assumes axis sizes are concomittant I think
    def toXarrayDS(self, N = None):

        if N is None:
            N = self.Emod[N] -1   # Default to last set field set

        ds = xr.Dataset()  # Init empty dataset, then loop over fields

        for domain in self.Emod[N].keys():

            domainName = f'E{domain}'  # Set data name - can't be same as coord names in this case

            if domain in [self.Edef['Pulse']['domain'], self.Edef['Pulse']['dfft']]:  # Key dims, set as unlinked
                ds[domainName] = ((domain), self.Emod[N][domain]['Ef'])
                ds.coords[domain] = self.Emod[N][domain]['axis']

            else:
                ds[domainName] = ((domain), self.Emod[N][domain]['Ef'])
                ds.coords[domain] = ((self.Edef['Pulse']['dfft']), self.Emod[N][domain]['axis'])  # For derived dims, set as linked to dfft dim - should set this more cleanly elsewhere...?


            # Set units (will be used for plotting)
            ds[domain].attrs['units'] = self.Edef['Units'][domain]['unit']
            ds[domainName].attrs['units'] = 'arb'

        # Assign spectrogram
        if hasattr(self, 'Spectrogram'):
            ds['spectrogram'] = ((self.Edef['Pulse']['dfft'], self.Edef['Pulse']['domain']), np.abs(self.Spectrogram['data']))


        self.Emod[N]['ds'] = ds



    #***************** Plotting

    # Basic plotting
    # TODO: check plotTypes vs. defns. in ePSproc
    # TODO: modify to plot sets of fields, either stacked or single plot.
    # 28/04/20: Changed to plot from Emod[N] dicts
    # TODO: sort out axis limits - should pass to override or set in class. Also change to, e.g. Holoviews, for more interaction...
    def plot(self, plotType = 'phaseUW', Nplot = None, domainList = None, thres = 1e-2):
        '''Basic plot with Matplotlib'''
#         plt.plot(self.t, self.Et)

        if Nplot is None:
            Nplot = np.arange(0, self.Emod['N'])

        # Default plots for domain + dfft
        if domainList is None:
            domainList = [self.Edef['Pulse']['domain'], self.Edef['Pulse']['dfft']]

        for domain in domainList:
            # Set up figure - do this *before* looping over N (fields)
            fig, ax1 = plt.subplots()
            if plotType == 'phaseUW':
                ax2 = ax1.twinx()

#             lText = []

            # Plot selected fields
            for N in Nplot:
                lString = f'$E_{{{N}}}$ '

                # Plot according to type
                if plotType == 'complex':
                    # Plot real + imag field components
                    ax1.plot(self.Emod[N][domain]['axis'], self.Emod[N][domain]['Ef'].real, '-', label = lString + 'Re')
                    ax1.plot(self.Emod[N][domain]['axis'], self.Emod[N][domain]['Ef'].imag, '-', label = lString + 'Im')
                    ax1.plot(self.Emod[N][domain]['axis'], np.abs(self.Emod[N][domain]['Ef']), '--', label = f'|{lString}|')
#                     plt.legend(['Re', 'Im', 'Abs'])
#                     lText.extend([f'{N} Re', f'{N} Im', f'{N} Abs'])

                elif plotType == 'field':
                    # Plot real-valued field (E + E*)
                    ax1.plot(self.Emod[N][domain]['axis'], 0.5*(self.Emod[N][domain]['Ef'] + self.Emod[N][domain]['Ef'].conj()), '-', label = f'{lString}+{lString}*')


                elif plotType == 'abs':
                    # Plot envelope only, |E|
                    ax1.plot(self.Emod[N][domain]['axis'], np.abs(self.Emod[N][domain]['Ef']), '-', label = f'|{lString}|')


                elif plotType == 'phase':
                    # Plot magnitude + phase
                    ax1.plot(self.Emod[N][domain]['axis'], np.abs(self.Emod[N][domain]['Ef']),'-', label = f'|{lString}|')
                    ax1.plot(self.Emod[N][domain]['axis'], (np.angle(self.Emod[N][domain]['Ef'])), '--', label = lString + 'Phase')
#                     lText.extend([f'{N} Abs', f'{N} Phase'])
#                     lText.extend((f'{N} Abs', f'{N} Phase'))
#                     ax1.legend(['Abs', 'Phase'])

                elif plotType == 'phaseUW':
                    # Plot magnitude + phase, unwrapped

                    # Single axis
    #                 plt.plot(self.Edef[domain]['axis'], np.abs(self.Edef[domain]['Ef']), self.Edef[domain]['axis'], np.unwrap(np.angle(self.Edef[domain]['Ef'])))

                    # Test secondary_y - not working
    #                 plt.plot(self.Edef[domain]['axis'], np.abs(self.Edef[domain]['Ef']))
    #                 plt.plot(self.Edef[domain]['axis'], np.unwrap(np.angle(self.Edef[domain]['Ef'])), secondary_y=True)

                    # Full ax addressing
                    ax1.plot(self.Emod[N][domain]['axis'], np.abs(self.Emod[N][domain]['Ef']), '-', label = f'|{lString}|')
                    ax2.plot(self.Emod[N][domain]['axis'], np.unwrap(np.angle(self.Emod[N][domain]['Ef'])), '--', label = lString + 'Phase')

#                     plt.legend(['Abs', 'Phase (unwrapped)'])
#                     lText.extend([f'{N} Abs', f'{N} Phase'])
#                     lText.extend((f'{N} Abs', f'{N} Phase'))


            if plotType != 'phaseUW':
                plt.ylabel('Amplitude')
                plt.xlabel(self.Edef['Units'][domain]['unit'])

            else:
                ax1.set_ylabel('Amplitude')
                ax2.set_ylabel('Phase (unwrapped)')
                ax1.set_xlabel(self.Edef['Units'][domain]['unit'])


        #         plt.xlim((-2.66, 2.66))  # Set for two cycles at 800nm
    #             plt.xlim(-0.5, 0.5)

            # Set some sensible limits (FFT scales will be large)
            if domain == 'f':  # self.Edef['Pulse']['dfft']:  # Hmmm, should be able to do this generically over all domains?  With origin + width?
                                                                # TODO: move origins + widths to domain-specific containers.
                plt.xlim(0.8*self.Edef['Freq']['f'], 1.2*self.Edef['Freq']['f'])

            elif domain == 't':
                scale = [-2.5, 2.5]
                plt.xlim(scale[0]*self.Edef['Pulse']['FWHM'], scale[1]*self.Edef['Pulse']['FWHM'])

            else:
                # Estimate from feature...
                peak = np.abs(self.Emod[N][domain]['Ef']).max()
                mask = np.abs(self.Emod[N][domain]['Ef']) > peak*thres
                scale = [self.Emod[N][domain]['axis'][mask].min(), self.Emod[N][domain]['axis'][mask].max()]
                plt.xlim(scale[0], scale[1])

            # Set legend from array or list
#             plt.legend(Nplot)
#             plt.legend(lText)

            # Set legends from labels, per axis object
            if plotType == 'phaseUW':
                ax1.legend(loc='upper left')
                ax2.legend(loc='upper right')
            else:
                ax1.legend(loc='upper left')  # This sometimes defaults to middle, so set explicitly


            plt.title(self.Estring + f'\nplotType = {plotType}')
            plt.show()


    def plotSpectrogram(self):
        """
        VERY basic spectrogram plotter from old code - for quick testing only.

        TODO: fix axes for with/withou fft shift.

        BETTER: use plotSpectrogramHV() instead, this uses full axes as set.

        """

        # Testing - set vars as per old code for brevity
        S = self.Spectrogram['data']
        N = self.Emod['N'] - 1
        t = self.Emod[N]['t']['axis']
        f = self.Emod[N]['f']['axis']

        plt.figure()
        plt.imshow(np.abs(S)**2, extent = [t[0],t[-1],f[np.int(S.shape[0]/2)-1],f[np.int(S.shape[0]/2)]], aspect='auto')
        plt.ylim((1.5*(f0/tUnit),2.5*(f0/tUnit)))
        plt.ylabel('$\Omega$ /THz')
        plt.xlabel('t /fs')
        plt.title('Spectrogram')
        plt.show()



    def plotSpectrogramHV(self, N = None):
        # May need to set this directly in notebook?
        hv.extension('bokeh')

        # Set N
        if N is None:
            N = self.Emod['N'] -1

        # Check ds exists
        if not 'ds' in self.Emod[N].keys():
            self.toXarrayDS(N = N)


        hv_ds = hv.Dataset(self.Emod[N]['ds']['spectrogram'])
#         hv_ds
        self.spec = hv_ds.to(hv.Image, kdims=['t','f'])
        self.spec.opts(width=700, height=700)
