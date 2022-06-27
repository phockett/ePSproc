"""
ePSproc io util functions

Various tools for use in Xarray file IO.

27/06/22    Split out from core IO.py, to extend backend options and support.
            Now additionally wrapped therein for flexible handling of multiple backends.

"""

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
