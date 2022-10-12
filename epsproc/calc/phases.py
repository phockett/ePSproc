"""
ePSproc phase calculations and wrappers

12/10/22    v1 basic implementation.

"""

# Imports
import numpy as np
import xarray as xr
import mpmath as mp   # Used for Gamma functions

from epsproc.util.conversion import conv_ev_atm

def coulombPhase(data = None, dims = {'E':'Eh','l':'l'},
                 lList = None, eList = None, eConv = False,
                 Z = [1, 1],
                 returnType = 'xr'):
    r"""
    Calculate Coulomb phases

    .. math::

       \sigma_{l,k}=\arg\Gamma\left[l+1-i\frac{Z_{1}Z_{2}}{k}\right]

    Where :math:`\Gamma` is the gamma function.

    For more details, see `Coulomb wavefunction notes <https://github.com/phockett/Quantum-Metrology-with-Photoelectrons/blob/master/Coulombic_case/Continuum_wavefunctions.ipynb>`_


    Parameters
    ----------

    data : named data array, e.g. Xarray or Pandas, default = None
        Data array containing coords l and E to use.
        Any datatype which supports dict retrieval, e.g. data['l'], data['E'] can be passed here.

    dims : dict, optional, default = {'E':'Eh','l':'l'}
        For named array case, this can be passed to override default dims (l,Eh).
        For Xarray outputs the same dim names will be used.
        E.g. dims = {'E':'Eke','l':'L'} will use passed data dims (Eke,L).

    lList : list or array, optional, default = None
        List of l values to calculate for.

    eList : list or array, optional, default = None
        List of E/k values to calcualte for.
        Note these must be in atomic units (Hartree).

    eConv : bool, default = False
        If true convert eList from eV to Hartrees.
        If Xarray data is passed with units set this will be forced to True if data[dims['E']].attrs['units'] == 'eV'

    Z : list, default = [1,1]
        Defines values [Z1,Z2].

    returnType : str, optional, default = 'xr'
        Return data as
        - 'xr' for Xarray
        - 'np' for Numpy array


    Note
    ----
    Pass EITHER data (plus optional dims) OR lList, eList values.

    Uses mpmath library for :math:`\Gamma`, see https://mpmath.org/doc/current/functions/gamma.html#gamma

    """

     # Set default dim names, also used for XR output.
#     if dims is None:
#         dims = {'E':'Eh','l':'l'}

    # Init from passed data array
    if data is not None:

        # TODO: add dim check here?

        lList = np.unique(data[dims['l']])  # Use unique here always?  Stacked LM index may have duplicates.
        eList = data[dims['E']]

        # For Xarray need to pull values.
        if hasattr(eList,'values'):
            eList = eList.values

        # For Xarray may have units listed, and convert if required
        if hasattr(data[dims['E']], 'attrs'):

            if 'units' in data[dims['E']].attrs.keys():
                if data[dims['E']].attrs['units'] == 'eV':
#                     eList = ep.conv_ev_atm(eList, to='H')
                    eConv = True

#                 else:
#                     eConv = False


    if eConv:
        eList = conv_ev_atm(eList, to='H')



    # Compute phases
    # For mpmath need to loop over all values & convert outputs.
    cPhase = {}
    sigmaNP = []

    for n,k in enumerate(eList):
        cPhase[n] = {'k':k,'sigmaMP':[],'l':[]}     # ,'l':lList}  # Set list in loop to ensure order(?)

        for l in lList:
            cPhase[n]['l'].append(l)
            cPhase[n]['sigmaMP'].append([mp.arg(mp.gamma(l+1-((1j*Z[0]*Z[1])/k)))])
    #         cPhase[n]['sigma'].append([l, mp.arg(mp.gamma(l+1-(1j/k)))])

        sigmaNP.append(cPhase[n]['sigmaMP'])

    # Stack outputs
    sigmaNP = np.asarray(sigmaNP,dtype=float).squeeze()   # float or np.float128?

    if returnType == 'np':
        return sigmaNP

    else:
        # TODO - add some attrs here.
        return xr.DataArray(sigmaNP, coords = {dims['E']:eList, dims['l']:lList})


def phaseCorrection(data, cPhase = None, lPhase = True, **kwargs):
    r"""
    Calculate and apply dipole phase correction term to a data array:

    .. math::
        i^{-l}*exp[i \sigma_l(k)]

    Parameters
    ----------

    data : Xarray
        Array to apply phase correction to

    cPhase : Xarray, optional, default = None
        Coulomb phases to use.
        If None, compute these for input array using :py:func:`coulombPhase()` function.
        **kwargs are also passed to :py:func:`coulombPhase()`

    lPhase : bool, default = True
        Apply :math:`i^{-l}` term?
        Skip this term if false

    **kwargs : optional keyword args
        Passed to :py:func:`coulombPhase()`
        

    Returns
    -------

    Xarray : data * phaseCorr

    Xarray : phaseCorr


    """

    if cPhase is None:
        cPhase = coulombPhase(data, **kwargs)

    if lPhase:
        lPhase = 1/(1j**cPhase.l)   # i^(-l) == 1/(i^l)

    else:
        lPhase = 1

    phaseCorr = lPhase*np.exp(1j*cPhase)

    # TODO: may need to check l dim and unstack here?
    return data * phaseCorr, phaseCorr
