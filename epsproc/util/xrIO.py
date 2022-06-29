"""
ePSproc Xarray IO util functions

Various tools for use in Xarray file IO.

27/06/22    Split out from core IO.py, to extend backend options and support.
            Now additionally wrapped therein for flexible handling of multiple backends.

"""

import numpy as np
import xarray as xr

#*********** Complex data handling

# Split complex to R + I
def splitComplex(data):
    """Split complex data into R+I floats."""

    dataR = np.real(data)
    dataI = np.imag(data)

    return dataR, dataI

# Comibine R + I to complex
def combineComplex(dataR, dataI):
    """Combine R+I floats into complex form."""

    data = dataR + 1j*dataI

    return data


def splitComplexXR(dataIn):
    """
    Split complex-valued Xarray data & coords to Re + Im components

    Splits input Xarray into Xarray Dataset with 'Re' and 'Im' components.

    """

    # Safe version with re/im split save type only.
    # Works for scipy and h5netcdf OK, latter will save complex type too, but is strictly not valid.
    dataOut = xr.Dataset({'Re':dataIn.real, 'Im':dataIn.imag})
    # dataOut.attrs = dataIn.attrs   # This will push dataarray attrs to dataset attrs, otherwise they're nested
                                    # May not always want this?

    # Allow for SF & XS coords which may also be complex
    # if 'XS' in dataOut.coords:
    #     dataOut['XSr'] = dataOut.XS.real
    #     dataOut['XSi'] = dataOut.XS.imag
    #     dataOut = dataOut.drop('XS')
    #
    # if 'SF' in dataOut.coords:
    #     dataOut['SFr'] = dataOut.SF.real
    #     dataOut['SFi'] = dataOut.SF.imag
    #     dataOut = dataOut.drop('SF')

    # Allow for arb complex coords.
    # May also want to add attr checker here? Or set in 'sanitizeAttrsNetCDF'
    for item in dataOut.coords.keys():
        if dataOut.coords[item].dtype == 'complex128':
            dataOut.coords[item + 'r'], dataOut.coords[item + 'i'] = splitComplex(dataOut.coords[item])
            dataOut = dataOut.drop(item)

    # Force top-level attrs
    dataOut.attrs = dataOut[list(dataOut.data_vars)[0]].attrs
    dataOut.attrs['complex'] = 'split'

    return dataOut


def combineComplexXR(dataIn):
    """
    Combine Re + Im Xarray Dataset and coordinates to complex values

    Note: not general, assumes formatting as defined by splitComplexXR()

    """

    # Reconstruct complex variables, NOTE this drops attrs... there's likely a better way to do this!
    # UPDATE 07/06/22: additional attrs handling below. Note in this case dataOut is a DataArray here.
    dataOut = dataIn.Re + dataIn.Im*1j
    # dataOut.attrs = dataIn.attrs

    # Rest SF & XS coords which may also be complex
    # Note: need to check vs. dataIn here, since dataOut already has dropped vars
    # if 'XSr' in dataIn.data_vars:
    #     dataOut['XS'] = dataIn.XSr + dataIn.XSi*1j
    # #     dataOut = dataOut.drop('XSr').drop('XSi')
    #
    # if 'SFr' in dataIn.data_vars:
    #     dataOut['SF'] = dataIn.SFr + dataIn.SFi
    # #     dataOut = dataOut.drop('SFr').drop('SFi')

    # General version
    for item in dataOut.coords.keys():
        # Check for r+i pairs - note labelling assumed to match writeXarray conventions here.
        if item.endswith('r'):
            itemi = item[:-1] + 'i'

            # If imag partner found, restack and remove split components.
            if itemi in dataOut.coords.keys():
                dataOut.coords[item[:-1]] = combineComplex(dataOut.coords[item], dataOut.coords[itemi])
                dataOut = dataOut.drop([item,itemi])

    # For dataset case, try some generic handling. May need more sophisticated methods here, maybe just assume DataArray and convert?
    if (not dataOut.attrs) and isinstance(dataIn, xr.core.dataset.Dataset):
        dataOut.attrs = dataIn[list(dataIn.data_vars)[0]].attrs
    # dataOut.attrs = dataIn[list(dataIn.data_vars)[0]].attrs

    return dataOut


#*********** Attribs handling

# Sanitize attributes & dicts for Xarray NetCDF IO
def sanitizeAttrsNetCDF(data, dictHandling = 'wrap'):
    """
    Sanitize Xarray DataArray attributes for file IO.

    Note this may be lossy:

    - Empty data > string.
    - Dictionaries removed, wrapped to string, or left alone (nested dicts not supported in attrs for most (all?) file writers).
      Set dictHandling = 'del', 'wrap' or anything else to leave as is.
    - Remove all items not of types [str, np.ndarray, int, float, list, tuple]


    Todo:

    - try conversion to string for all attrs?
    - try dict conversions & JSON side-car file IO to avoid lossy saves.

    """

    dataOut = data.copy()

    # Remove None and other empty types, ugh - now integrated below
    # xrTest.attrs = {k:(v if v else str(v)) for k,v in xrTest.attrs.items()}
    log = {}
    for k,v in dataOut.attrs.items():
        if not v:
            dataOut.attrs[k] = str(v)
            log[k] = 'str'

        if isinstance(dataOut.attrs[k], dict):
    #         xrTest.attrs[k] = [[k2,v2] for k2,v2 in xrTest.attrs[k].items()]  # Nest dict items also not supported, dump to nested lists? Seems to be acceptable. Ugh.
                                                                                # Still causing issues in some cases?
            if dictHandling == 'del':
                dataOut.attrs[k] = 'Removed dict'
                log[k] = 'Removed dict'

            elif dictHandling == 'wrap':
                dataOut.attrs[k] = str(v)
                log[k] = 'Wrapped dict to string'

            else:
                pass

        if type(dataOut.attrs[k]) not in [str, np.ndarray, int, float, list, tuple]:
            typeIn = type(dataOut.attrs[k])
            dataOut.attrs[k] = 'NA'
            log[k] = f'Removed item type {typeIn}'

    # TO TRY - full str conversion, e.g. from https://stackoverflow.com/a/42676094 (for JSON example case)
    # save: convert each tuple key to a string before saving as json object
    # s = json.dumps({str(k): str(v) for k, v in eulerDict.items()})
    #
    # THEN RECON with ast:
    # # load in two stages:
    # # (i) load json object
    # obj = json.loads(s)
    #
    # # (ii) convert loaded keys from string back to tuple
    # from ast import literal_eval
    # # d = {literal_eval(k): literal_eval(v) for k, v in obj.items()}  # FAILS: ValueError: malformed node or string: <ast.Name object at 0x7f4464e67550>
    # d = {k: (literal_eval(v) if v != 'Euler' else v) for k, v in obj.items()}  # ok - WORKS FOR ALL CASES EXCEPT NON-EXECUTABLE STRS
    #
    # This should also work here, but maybe add type checking too?


    return dataOut, data.attrs, log
