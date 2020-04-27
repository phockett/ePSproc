"""
ePSproc electric field functions.

Generate, propagate and model E-fields.


27/04/20    Updating & packaging - but still very much a work-in-progress.

20/3/20     Basics developed in tests/methodDev/geometric_method_dev_low-level_E-fields_200320.ipynb


"""


import numpy as np
from scipy import constants as scipy_constants
import pprint

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

    '''

    def __init__(self, Edef = None, units = None):  # , E0 = None, sigma = None, A = 1, t = None, dt = None, l0 = None, f0 = None, w0 = None, units = None):

        #*** Define master list of parameters & init as empty or with defaults
        # If no parameters are passed, this can be used to create a defn template
        self.Edef = {'Pulse' : {'type' : 'Gaussian',   # Eventually should generalise this to allow for other pulse types.
                                'domain' : 't',       # Set a domain label here to switch between (t,E) domain pulse defn. - MAY WANT TO MOVE THIS TO SEPARATE LIST FN?
                                'dfft' : 'f',        # Set conjugate domain (for FFT)
                                'sigma' : None,
                                'E0' : 1,
                                'A' : 1,
                                'FWHM' : None,
                                'p' : 0    # Set polarization state. For single values {-1,0,+1} this will just label output. For 'XY' calc. spherical terms.
#                                 'origin' : 0    # Set origin for domain
                               },

                     'Freq' : {'l' : None,   # Set defns. for carrier, (wavelength, freq, ang. freq). Assumed to be in working units. Also use to index units later.
                               'f' : None,
                               'w' : None,
                               # 'Ehv' : None  # Should add here, but need to fix assignment routine first!
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

                     'Spectrogram' : {'gate' : None,
                                      'data' : None
                                     },

                     'Units' : {   # 'f':{'value':1e12, 'unit':'THz'},  # Set this as a derived unit, 1/t, below, otherwise may get discrepancies in unit conversion.
                                'l':{'value':1e-9, 'unit':'nm'},
                                't':{'value':1e-15, 'unit':'fs'},
                                'Ehv':{'value':1, 'unit':'eV'}
                                },

                     # Global FFT settings
                     'FFT' : {'pad':True, 'positiveHalf':True, 'phaseMaskFlag':False, 'thres':1e-3}
                    }

        self.Emod = {}  # Set empty dict to hold modified fields (propagated, shaped etc.)


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

            # Set units for E
            self.Edef['Units']['E'] = self.Edef['Units']['f'].copy()

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
                    self.Edef['t']['delta'] = 1e-3*self.Edef['Pulse']['sigma']  # Step size for t axis, relative to sigma. May get large however!

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
                        print(refValue)

                if refKey is 'w':
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

                if refKey is not 'w':  # Set w0 in working units (rad/[unit t])
#                     self.Edef['Freq']['w'] = 2*np.pi * self.Edef['Freq']['f'] * self.Edef['Units']['f']['value'] * self.Edef['Units']['t']['value']
                    self.Edef['Freq']['w'] = 2*np.pi * self.Edef['Freq']['f']


            self.ECalc()   # Set defined field
            self.EFFT()    # Set FFT field


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
            self.Edef[domain]['Ef'] = self.Edef[domain]['Ef'] * np.exp(1.0j*self.Edef['Freq']['w']*self.Edef[domain]['axis'])


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

        Ebar = np.fft.fft(self.Edef[domain]['Ef'], n=n)
        axis = np.fft.fftfreq(nAxis, d=self.Edef[domain]['delta'])

        # Set for full or half FT
        if self.Edef['FFT']['positiveHalf']:
            inds = np.int(Ebar.shape[0]/2)   # Set for half the range
            self.Edef[dfft]['Ef'] = Ebar[0:inds]
            self.Edef[dfft]['axis'] = axis[0:inds]
            self.Edef[dfft]['delta'] = axis[1]-axis[0]
        else:
            self.Edef[dfft]['Ef'] = Ebar
            self.Edef[dfft]['axis'] = axis
            self.Edef[dfft]['delta'] = axis[1]-axis[0]

        if self.Edef['FFT']['phaseMaskFlag']:
            self.phaseMask(domain = dfft)


    # Calculate iFFT(E).
    # This, sort of, assumes that the field is defined as E(w), and the spectral phase is modified.
    # But could be used for other cases.
    # Emod is used if passed, otherwise uses self.Edef[dfft]['Ef']
    # 27/04/20 Modified to use self.Emod for fields.
    def EiFFT(self, Emod = None): #, f=None, pad=False):

        domain = self.Edef['Pulse']['domain']
        dfft = self.Edef['Pulse']['dfft']

        if Emod is None:
            Emod = self.Edef[dfft]['Ef']

        # Transform back to time-domain, and center
        Eifft = np.fft.ifftshift(np.fft.ifft(Emod))   # With shift
#         Eifft = (np.fft.ifft(self.Edef[dfft]['Ef']))  # Without shift

        if self.Edef['FFT']['positiveHalf']:  # Correct amplitudes, get x2 otherwise if positiveHalf is True.
            Eifft /= 2.0



#         self.Edef[domain]['Ef'] = np.c_[self.Edef[domain]['Ef'], Eifft]   # Basic stacking OK if pad=False and positiveHalf=False, otherwise axis lengths different

        # Set for full or half FT - if used with iFFT shift and a pulse center at 0, this will slice result.
        inds = np.int(Eifft.shape[0])   # Set for full range
#         if self.Edef['FFT']['positiveHalf']:
#             inds = np.int(Eifft.shape[0]/2)   # Set for half the range

        # Set as new field, or as new dict?
        # Should change formatting here, may want domain_N to allow for multiple cases, and looping.
        if (domain + '_mod') in self.Edef.keys():
            N = list(self.Edef[domain + '_mod'].keys())
        else:
            self.Edef[domain + '_mod'] = {}
            N = 0

        self.Edef[domain + '_mod'][N] = {'Ef':Eifft[0:inds],
                                     'axis':np.fft.ifftshift(np.fft.fftfreq(Eifft.shape[-1], d=self.Edef[dfft]['delta'])[0:inds]),
                                     'domain':domain

                                    }


        #if f:
        #    t = np.fft.ifftfreq(f.shape[-1])
        #else:
        #    t = None


    #*************** Phase modification fns.

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



    def chirp(self, A, resetPhase=True):
        """
        Add quadratic phase (chirp) in spectral domain. Requires an E-field object, and chirp parameter A.

        .. :math: phi=A*(f-f0)^2

        TODO: check Diels for nomenclature here.

        """
        domain = 'f'

        # Remove existing phase - this is only the linear term if starting from a Gaussian pulse
        if resetPhase:
            Ew = np.abs(self.Edef[domain]['Ef'])
        else:
            Ew = self.Edef[domain]['Ef']

        # Add chirp
#         Ewq = np.abs(Ew)*np.exp(1.0j*(A*(self.Edef[domain]['axis']-self.Edef['Freq']['f']))**2)
        Ewq = np.abs(Ew)*np.exp(1.0j*A*(self.Edef[domain]['axis']-self.Edef['Freq']['f'])**2)


        # Set as new field, or as new dict?
        # Should change formatting here, may want domain_N to allow for multiple cases, and looping.
        # Or nest this? Or concatenate? (-- Only if domain axes are identical.)
        if (domain + '_mod') in self.Edef.keys():
            N = list(self.Edef[domain + '_mod'].keys())
        else:
            self.Edef[domain + '_mod'] = {}
            N = 0

        self.Edef[domain + '_mod'][N] = {'Ef':Ewq,
                                     'axis':self.Edef[domain]['axis'],
                                     'domain':domain

                                    }
#         self.Edef[domain] = {'Ef':Ewq,
#                                      'axis':self.Edef[domain]['axis'],
#                                      'domain':domain

#                                     }

        self.EiFFT(Emod = Ewq)  # Set tmod pulse via ifft



    #***************** Plotting

    # Basic plotting
    # TODO: check plotTypes vs. defns. in ePSproc
    def plot(self, plotType = 'phaseUW'):
        '''Basic plot with Matplotlib'''
#         plt.plot(self.t, self.Et)

        for domain in [self.Edef['Pulse']['domain'], self.Edef['Pulse']['dfft']]:

            # Plot according to type
            if plotType == 'complex':
                # Plot real + imag field components
                plt.plot(self.Edef[domain]['axis'], self.Edef[domain]['Ef'].real, self.Edef[domain]['axis'], self.Edef[domain]['Ef'].imag, self.Edef[domain]['axis'], np.abs(self.Edef[domain]['Ef']))
                plt.legend(['Re', 'Im', 'Abs'])

            elif plotType == 'phase':
                # Plot magnitude + phase
                plt.plot(self.Edef[domain]['axis'], np.abs(self.Edef[domain]['Ef']), self.Edef[domain]['axis'], (np.angle(self.Edef[domain]['Ef'])))
                plt.legend(['Abs', 'Phase'])

            elif plotType == 'phaseUW':
                # Plot magnitude + phase, unwrapped

                # Single axis
#                 plt.plot(self.Edef[domain]['axis'], np.abs(self.Edef[domain]['Ef']), self.Edef[domain]['axis'], np.unwrap(np.angle(self.Edef[domain]['Ef'])))

                # Test secondary_y - not working
#                 plt.plot(self.Edef[domain]['axis'], np.abs(self.Edef[domain]['Ef']))
#                 plt.plot(self.Edef[domain]['axis'], np.unwrap(np.angle(self.Edef[domain]['Ef'])), secondary_y=True)

                # Full ax addressing
                fig, ax1 = plt.subplots()
                ax2 = ax1.twinx()
                ax1.plot(self.Edef[domain]['axis'], np.abs(self.Edef[domain]['Ef']), 'g-')
                ax2.plot(self.Edef[domain]['axis'], np.unwrap(np.angle(self.Edef[domain]['Ef'])), 'b-')

                plt.legend(['Abs', 'Phase (unwrapped)'])

            if plotType != 'phaseUW':
                plt.ylabel('Amplitude')
                plt.xlabel(self.Edef['Units'][domain]['unit'])

            else:
                ax1.set_ylabel('Amplitude')
                ax2.set_ylabel('Phase (unwrapped)')
                ax1.set_xlabel(self.Edef['Units'][domain]['unit'])


    #         plt.xlim((-2.66, 2.66))  # Set for two cycles at 800nm
#             plt.xlim(-0.5, 0.5)

            if domain == 'f':  # self.Edef['Pulse']['dfft']:  # Hmmm, should be able to do this generically over all domains?  With origin + width?
                                                                # TODO: move origins + widths to domain-specific containers.
                plt.xlim(0.8*self.Edef['Freq']['f'], 1.2*self.Edef['Freq']['f'])

            plt.title(self.Estring)
            plt.show()
