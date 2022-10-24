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
def setDataPath(setDemoDataPath):
    # Set data path
    # Note this is set here from ep.__path__, but may not be correct in all cases - depends on where the Github repo is.
    # epDemoDataPath = Path(ep.__path__[0]).parent/'data'
    epDemoDataPath = setDemoDataPath

    # dataPath = os.path.join(epDemoDataPath, 'photoionization', 'n2_multiorb')
    dataPath = Path(epDemoDataPath, 'photoionization', 'n2_multiorb')

    return dataPath


@pytest.fixture(scope="module")  # Scoping here doesn't help with class issues below.
def setAlignmentDataFile(setDemoDataPath):
    # Set data path
    # Note this is set here from ep.__path__, but may not be correct in all cases - depends on where the Github repo is.
    # epDemoDataPath = Path(ep.__path__[0]).parent/'data'
    epDemoDataPath = setDemoDataPath

    # dataPath = os.path.join(epDemoDataPath, 'photoionization', 'n2_multiorb')
    dataPath = Path(epDemoDataPath, 'alignment', 'N2_ADM_VM_290816.mat')

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


# Test AFBLM calcs for aligned case.
# Adapted from https://pemtk.readthedocs.io/en/latest/fitting/PEMtk_fitting_basic_demo_030621-full_010922.html
def test_AFBLM_N2_ADMs(dataClass, setAlignmentDataFile):

    data = dataClass

    #*** Load time-dependent ADMs for N2 case
    # Adapted from ePSproc_AFBLM_testing_010519_300719.m

    from scipy.io import loadmat
    # ADMdataFile = os.path.join(epDemoDataPath, 'alignment', 'N2_ADM_VM_290816.mat')
    ADMdataFile = setAlignmentDataFile
    ADMs = loadmat(ADMdataFile)

    # Set tOffset for calcs, 3.76ps!!!
    # This is because this is 2-pulse case, and will set t=0 to 2nd pulse (and matches defn. in N2 experimental paper)
    tOffset = -3.76
    ADMs['time'] = ADMs['time'] + tOffset


    # BASE VERSION FROM https://epsproc.readthedocs.io/en/dev/methods/geometric_method_dev_pt3_AFBLM_090620_010920_dev_bk100920.html#Test-compared-to-experimental-N2-AF-results...
    # Run with sym summation...
    phaseConvention = 'E'  # Set phase conventions used in the numerics - for ePolyScat matrix elements, set to 'E', to match defns. above.
    symSum = True  # Sum over symmetry groups, or keep separate?
    SFflag = False  # Include scaling factor to Mb in calculation?

    SFflagRenorm = False  # Renorm terms
    BLMRenorm = 1

    thres = 1e-4
    # ADMX = ep.setADMs(ADMs = ADMs['ADM'], t=ADMs['time'].squeeze(), KQSLabels = ADMs['ADMlist'], addS = True)
    RX = ep.setPolGeoms()  # Set default pol geoms (z,x,y), or will be set by mfblmXprod() defaults - FOR AF case this is only used to set 'z' geom for unity wigner D's - should rationalise this!


    # Selection & downsampling
    trange=[4, 5]  # Set range in ps for calc
    tStep=4  # Set tStep for downsampling

    tMask = (ADMs['time']>trange[0]) * (ADMs['time']<trange[1])
    ind = np.nonzero(tMask)[1][0::tStep]
    At = ADMs['time'][:,ind].squeeze()
    ADMin = ADMs['ADM'][:,ind]

    ADMX = ep.setADMs(ADMs = ADMs['ADM'][:,ind], t=At, KQSLabels = ADMs['ADMlist'], addS = True)

    # Original version with base function
    # start = time.time()
    # mTermST, mTermS, mTermTest, BetaNormX = ep.geomFunc.afblmXprod(dataSet[0], QNs = None, AKQS=ADMX, RX=RX, thres = thres, selDims = {'it':1, 'Type':'L'}, thresDims='Eke',
    #                                                     symSum=symSum, SFflag=SFflag, BLMRenorm = BLMRenorm,
    #                                                     phaseConvention=phaseConvention)
    # end = time.time()
    # print('Elapsed time = {0} seconds, for {1} energy points, {2} polarizations, threshold={3}.'.format((end-start), mTermST.Eke.size, RX.size, thres))

    # Class version.
    # NOte this currently has selection issues - dim needs to be (and remain) in matE, so currently fails for singleton Eke case
    start = time.time()
    data.AFBLM(AKQS=ADMX, RX=RX, thres = thres, selDims = {'it':1, 'Type':'L', 'Eke':[1.1,2.1]}, thresDims='Eke',   # thresDims='time',     # thresDims='Eke',
                                                    symSum=symSum, SFflag=SFflag, BLMRenorm = BLMRenorm,
                                                    phaseConvention=phaseConvention)
    end = time.time()

    # print('Elapsed time = {0} seconds, for {1} energy points, {2} polarizations, threshold={3}.'.format((end-start), mTermST.Eke.size, RX.size, thres))
    print(f'AFBLM calculations with alignment test, elapsed time = {(end-start)} seconds.')

    #*** Set data subsets
    # NOTE SOME OF THIS ONLY IMPLEMENTED IN PEMTK VERSION!
    # data.setADMs(ADMs = ADMs['ADM'], t=ADMs['time'].squeeze(), KQSLabels = ADMs['ADMlist'], addS = True)
    #
    # # Settings for type subselection are in selOpts[dataType]
    #
    # # Matrix element sub-selection
    # data.selOpts['matE'] = {'thres': 0.01, 'inds': {'Type':'L', 'Eke':1.1}}
    # data.setSubset(dataKey = 'orb5', dataType = 'matE')  # Subselect from 'orb5' dataset, matrix elements
    #
    # # Pols
    # # And for the polarisation geometries...
    # data.setPolGeoms()
    # data.selOpts['pol'] = {'inds': {'Labels': 'z'}}
    # data.setSubset(dataKey = 'pol', dataType = 'pol')
    #
    # # And for the ADMs...
    # data.selOpts['ADM'] = {}   #{'thres': 0.01, 'inds': {'Type':'L', 'Eke':1.1}}
    # data.setSubset(dataKey = 'ADM', dataType = 'ADM', sliceParams = {'t':[4, 5, 4]})
    #
    # #*** Compute
    # BetaNormX, basis = data.afblmMatEfit()  # OK, uses default polarizations & ADMs as set in data['subset']

    #*** SET DATA SUBSETS - basic version
    # See https://epsproc.readthedocs.io/en/dev/methods/geometric_method_dev_pt3_AFBLM_090620_010920_dev_bk100920.html#Test-compared-to-experimental-N2-AF-results...



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
