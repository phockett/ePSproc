# import os
# dataPath = os.path.join(sys.path[-1], 'data', 'photoionization', 'n2_multiorb')

import pytest
from pathlib import Path

import time
import numpy as np

# ePSproc
import epsproc as ep
from epsproc.classes.base import ePSbase

# @pytest.fixture(scope="session")
@pytest.fixture(scope="module")  # Scoping here doesn't help with class issues below.
def setDataPath():
    # Set data path
    # Note this is set here from ep.__path__, but may not be correct in all cases - depends on where the Github repo is.
    epDemoDataPath = Path(ep.__path__[0]).parent/'data'

    # dataPath = os.path.join(epDemoDataPath, 'photoionization', 'n2_multiorb')
    dataPath = Path(epDemoDataPath, 'photoionization', 'n2_multiorb')

    return dataPath


def test_base_init(setDataPath):

    data = ePSbase(setDataPath, verbose = 1)
    # data.scanFiles()

# Sequential load works
def test_base_scanFiles(setDataPath):

    data = ePSbase(setDataPath, verbose = 1)
    data.scanFiles()


#*** SCRATCH SPACE - testing methods for class setup and test

# Data as fixutre?
@pytest.fixture(scope="module")
def dataClass(setDataPath):

    data = ePSbase(setDataPath, verbose = 1)

    # Also load data
    data.scanFiles()

    return data

# Sequential load works
# def test_base_scanFiles2(dataClass):
#
#     # data = ePSbase(setDataPath, verbose = 1)
#     # data.scanFiles()
#     dataClass.scanFiles()


# Set return to match manually...
def jobsRef():
    refPrint = """
*** Job orb6 details
Key: orb6
Dir /home/femtolab/github/ePSproc/data/photoionization/n2_multiorb, 1 file(s).
{   'batch': 'ePS n2, batch n2_1pu_0.1-50.1eV, orbital A2',
    'event': ' N2 A-state (1piu-1)',
    'orbE': -17.096913836366,
    'orbLabel': '1piu-1'}

*** Job orb5 details
Key: orb5
Dir /home/femtolab/github/ePSproc/data/photoionization/n2_multiorb, 1 file(s).
{   'batch': 'ePS n2, batch n2_3sg_0.1-50.1eV, orbital A2',
    'event': ' N2 X-state (3sg-1)',
    'orbE': -17.341816310545997,
    'orbLabel': '3sg-1'}
"""

    return refPrint

# Test another method...
def test_jobsSummary(dataClass, capsys):
    # print('Testing jobs...')
    # print(dataClass.jobsSummary())
    dataClass.jobsSummary()

    captured = capsys.readouterr()
    assert captured.out == jobsRef()


# Test AFBLM calcultion for isotropic case
def test_AFBLM_iso(dataClass):

    data = dataClass

    # Calculate BLMs for default case
    start = time.time()
    data.AFBLM()
    end = time.time()
    print(f'AFBLM reference calc completed, elapsed time = {end-start} seconds')


    # Compare with ref ePS values to test
    # See https://epsproc.readthedocs.io/en/dev/methods/LF_AF_verification_tests_060720_tidy_100920.html#Test-vs.-AF-code

    # Note default threshold for AFBLM calcs, thres = 1e-2 - should give rough size of errors
    # May want to test for smaller thresholds too?
    thres = 1e-2

    # Set getCro options
    keys = data._keysCheck(None)
    pGauge = 'L'
    pTypes = [['SIGMA',0],['BETA',2]]  # Crude selection for XS and l values for test case only.
    sym = ('All','All')

    for key in keys:
        for pType in pTypes:

            if pType[0] == 'SIGMA':
                # Test XS
                XSdiff = data.data[key]['XS'].sel(XC=pType[0], Type=pGauge, Sym = sym) - data.data[key]['AFBLM'].sel(l=pType[1]).squeeze().XSrescaled  #, Type=pGauge)

            else:
                # Test B2, include conversion to Ylm coeffs.
                betaConv = ep.conversion.conv_BL_BLM(data.data[key]['XS'], to='sph', renorm=True)
                XSdiff = betaConv.sel(XC=pType[0], Type=pGauge, Sym = sym) - data.data[key]['AFBLM'].sel(l=pType[1]).squeeze()  #, Type=pGauge)


            assert np.abs(XSdiff.sum()) < (XSdiff.size * thres)
            assert np.abs(XSdiff.max()) < thres
            assert np.abs(XSdiff.max().imag)  < np.finfo(complex).eps   # Check complex zero to machine precision.


#******** IGNORE CLASS TESTS!!

# 2nd attempt at class test...
# Use setup() as per https://stackoverflow.com/a/39395889
# This seems to work... NO, same issues as below, fixture name not recognised?
# But trying to RUN fixture as `setDataPath()` yields error that fixture shouldn't be run directly.
# See also https://doc.pytest.org/en/latest/how-to/fixtures.html#usefixtures
# And https://docs.pytest.org/en/7.1.x/reference/fixtures.html#reference-fixtures

# pytest.mark.usefixtures(setDataPath)
# class TestBaseMethods():
#     def setup(self, setDataPath):
#         self.data = ePSbase(setDataPath, verbose = 1)
#          # self.data = ePSbase(setDataPath(), verbose = 1)
#
#     def test_scanFiles(self):
#         self.data.scanFiles()

# Class with global or local fixture fails?
# Getting incorrectly passed fixture somehow
#
# # # Best to test class methods with a class & data preloaded?
# class TestBaseMethods:
#     # data = ePSbase(setDataPath, verbose = 1)  # This fails? Fixture out of scope?
#
#     # def __init__():
#     #     data = ePSbase(setDataPath, verbose = 1)  # This fails, no __init__ allowed!
#
#     @pytest.fixture
#     def initClass(self):
#         self.data = ePSbase(setDataPath, verbose = 1)  # ALSO FAILS UPON RUN????
#         # return data
#
#
#     def test_scanFiles(self, initClass):
#
#         self.data.scanFiles()
