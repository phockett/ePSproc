"""
ePSproc verification functions

04/11/22    Added verifyAFiso(), based on test_AFBLM_iso, https://github.com/phockett/ePSproc/blob/4533a01c090a50ef69311e9718d32c3ce3dac9e7/tests/test_base_class.py#L106

"""
import time
import numpy as np

from epsproc.util.conversion import conv_BL_BLM

def verifyAFiso(data, thres = 1e-2, sym = ('All','All'), pGauge = 'L',
                outputDict = 'AFisoChecks'):
    """
    Verify AFBLM calculation for isotropic by comparison with ePS GetCro outputs.

    Based on test_AFBLM_iso() routine, https://github.com/phockett/ePSproc/blob/4533a01c090a50ef69311e9718d32c3ce3dac9e7/tests/test_base_class.py#L106

    Note that existing test routine for use with pyTest, this version for general use.

    Parameters
    ----------
    data : ePSproc data class object.

    thres : float, default = 1e-2

    sym : tuple, optional, default = ('All','All')
        Set symmetry selections for GetCro output.

    pGauge : str, optional, default = 'L'
        Set gauge selection fro GetCro output.

    outputDict : str, optional, default = 'AFisoChecks'
        Name for output data.
        Will be set to data.data[key][outputDict]

    Returns
    -------


    TODO: clean up input params, should match options elsewhere (e.g. matEleSelector), and use this for selections?

    """

    # data = dataClass

    # Calculate BLMs for default case
    start = time.time()
    data.AFBLM()
    end = time.time()
    print(f'\nAFBLM reference calc completed, elapsed time = {end-start} seconds')


    # Compare with ref ePS values to test
    # See https://epsproc.readthedocs.io/en/dev/methods/LF_AF_verification_tests_060720_tidy_100920.html#Test-vs.-AF-code

    # Note default threshold for AFBLM calcs, thres = 1e-2 - should give rough size of errors
    # May want to test for smaller thresholds too?
    # thres = 1e-2

    # Set getCro options
    keys = data._keysCheck(None)
    # pGauge = 'L'
    pTypes = [['SIGMA',0],['BETA',2]]  # Crude selection for XS and l values for test case only.
    # sym = ('All','All')

    outputSummary = {}

    for key in keys:

        data.data[key][outputDict] = {}
        outputSummary[key] = {}

        for pType in pTypes:

            if pType[0] == 'SIGMA':
                # Test XS
                XSdiff = data.data[key]['XS'].sel(XC=pType[0], Type=pGauge, Sym = sym) - data.data[key]['AFBLM'].sel(l=pType[1]).squeeze().XSrescaled  #, Type=pGauge)

            else:
                # Test B2, include conversion to Ylm coeffs.
                betaConv = conv_BL_BLM(data.data[key]['XS'], to='sph', renorm=True)
                XSdiff = betaConv.sel(XC=pType[0], Type=pGauge, Sym = sym) - data.data[key]['AFBLM'].sel(l=pType[1]).squeeze()  #, Type=pGauge)


            # assert np.abs(XSdiff.sum()) < (XSdiff.size * thres)
            # assert np.abs(XSdiff.max()) < thres
            # # assert np.abs(XSdiff.max().imag)  < np.finfo(complex).eps   # Check complex zero to machine precision.
            # assert np.abs(XSdiff.imag.max())  < np.finfo(complex).eps   # Check complex zero to machine precision.

            # Tests & outputs
            # Note .values or bool() to avoid full Xarray returns here.
            XSdiff.attrs['metrics'] = {'diffSum': np.abs(XSdiff.sum()).values,
                                       'diffSumOK': bool(np.abs(XSdiff.sum()) < (XSdiff.size * thres)),
                                       'diffAbsMax': np.abs(XSdiff.max()).values,
                                       'diffAbsMaxOK': bool(np.abs(XSdiff.max()) < thres),
                                       'diffImagMax': np.abs(XSdiff.imag.max()).values,
                                       'diffImagMaxOK': bool(np.abs(XSdiff.imag.max()) < np.finfo(complex).eps),
                                       }

            XSdiff.attrs['pass'] = XSdiff.attrs['metrics']['diffSumOK'] and XSdiff.attrs['metrics']['diffAbsMaxOK'] and XSdiff.attrs['metrics']['diffImagMaxOK']

            outputSummary[key]['metrics'] = XSdiff.attrs['metrics']
            outputSummary[key]['pass'] = XSdiff.attrs['pass']

            if XSdiff.attrs['pass']:
                print(f"AFBLM iso test for key={key}, dataType={pType} \t\t PASSED")
            else:
                print(f"AFBLM iso test for key={key}, dataType={pType} \t\t *** FAILED ***")

            data.data[key][outputDict][pType[0]] = XSdiff

    # Overall results
    outputSummary['All'] = all([outputSummary[key]['pass'] for key in outputSummary.keys()])

    if outputSummary['All']:
        print("\n*** All tests passed. ***")
    else:
        print("\n*** Some tests FAILED. ***")

    return outputSummary
