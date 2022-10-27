# Basic IO test routines.
# Mainly from https://epsproc.readthedocs.io/en/dev/dataStructures/ePSproc_dataStructures_IO_demo_280622.html

import pytest
from pathlib import Path
import pickle

# ePSproc
import epsproc as ep

# NOTE - should centralise this routine!
# UPDATE 23/10/22 - now moved fixtures to conftest.py

# # @pytest.fixture(scope="module")  # Scoping here doesn't help with class issues below.
# @pytest.fixture(scope="session")
# def setDataPath():
#     # Set data path
#     # Note this is set here from ep.__path__, but may not be correct in all cases - depends on where the Github repo is.
#     epDemoDataPath = Path(ep.__path__[0]).parent/'data'
#
#     # dataPath = os.path.join(epDemoDataPath, 'photoionization', 'n2_multiorb')
#     dataPath = Path(epDemoDataPath, 'photoionization')
#
#     return dataPath
#
# # @pytest.fixture(scope="module")
# @pytest.fixture(scope="session")
# def setDataFileName():
#     return 'n2_3sg_0.1-50.1eV_A2'
#
# # @pytest.fixture(scope="module")
# @pytest.fixture(scope="session")
# def setDataFile(setDataPath,setDataFileName):
#     dataPath = setDataPath
#     return Path(dataPath,setDataFileName).as_posix()  # Return POSIX to handle multiple suffixes later.
#                                                              # Or use .with_name?
#
# # Use tmp_path fixture for test data files, see https://docs.pytest.org/en/7.1.x/how-to/tmp_path.html
# # tmp_path_factory should be better for persistence? https://docs.pytest.org/en/7.1.x/how-to/tmp_path.html#the-tmp-path-factory-fixture
#
# # @pytest.fixture()
# # def setTmpDataFile(tmp_path,setDataFileName,tmp_path_factory):
# @pytest.fixture(scope="session")
# def setTmpDataFile(setDataFileName,tmp_path_factory):
#     # dataPath = tmp_path  # THIS WORKS, but note that this also sets tmp dirs PER TEST FUNCTION.
#     dataPath = tmp_path_factory.mktemp("data") / setDataFileName  # SAME ISSUE - sets DATAN per N run.
#                                                                   # NEEDS scope = "session" to use single dir over all tests.
#                                                                   # Note this may fail if other fixtures are NOT session scoped.
#     return dataPath.as_posix()
#
#     # return Path(dataPath,setDataFileName).as_posix()  # Return POSIX to handle multiple suffixes later.
#                                                              # Or use .with_name?
#
# # @pytest.fixture(scope="session")
# # def image_file(tmp_path_factory):
# #     img = compute_expensive_image()
# #     fn = tmp_path_factory.mktemp("data") / "img.png"
# #     img.save(fn)
# #     return fn
#
#
# # Data as fixutre?
# # @pytest.fixture(scope="module")
# @pytest.fixture(scope="session")
# # def dataSingle(setDataPath):
# def dataSingle(setDataFile):
#     # dataPath = setDataPath
#     # dataFile = Path(dataPath, 'n2_3sg_0.1-50.1eV_A2.inp.out')  # Set for sample N2 data for testing
#
#     # dataFile = setDataFile.with_suffix('.inp.out')
#     dataFile = setDataFile+'.inp.out'
#
#     # Scan data file
#     # dataSet = ep.readMatEle(fileIn = dataFile.as_posix())
#     dataSet = ep.readMatEle(fileIn = dataFile)
#     data = dataSet[0]
#
#     return data


def test_ePS_read(setDataPath, capsys, dataSingle, setDataFile):
    # Load data from modPath\data
    # dataPath = setDataPath
    # dataFile = Path(dataPath, 'n2_3sg_0.1-50.1eV_A2.inp.out')  # Set for sample N2 data for testing
    dataFile = setDataFile  # +'.inp.out'  # One liner OK in fixture, but not in function - execution order issues.
    dataFile = dataFile +'.inp.out'

    # Scan data file
    # dataSet = ep.readMatEle(fileIn = dataFile.as_posix())
    dataSet = ep.readMatEle(fileIn = dataFile)
    data = dataSet[0]

    # Fixture as test - should work, but need to return capsys output!
    # data = dataSingle

    captured = capsys.readouterr()
    # assert captured.out == loadRef()   # Can't get this to work for basic string case, issues with loadRef() formatting.
    # print(type(captured.out))
    # print(captured.out)

    # Write ref case to file
    # with open("ePS_read_ref.txt", "w") as text_file:
    #     text_file.write(captured.out)

    # Test vs. file - this works, now implemented as loadRef()
    # with open("Output.txt", "r") as text_file:
    # #     text_file.write(captured.out)
    #     ref = text_file.read()
    #
    # assert captured.out == ref

    # assert captured.out == loadRef(fileName = 'tests/ePS_read_ref.txt')
    # assert captured.out == loadRef(fileName = 'tests/ePS_read_ref_GHCI.txt')  # With paths set for GH action CI. TODO: don't test paths!
    # print('OK')


def test_pickle_write(dataSingle, setTmpDataFile):

    data = dataSingle

    # Save a dictionary
    # with open(Path(dataPath, 'n2_3sg_0.1-50.1eV_A2_dict.pickle'), 'wb') as handle:
    #     pickle.dump(dataDict, handle, protocol=pickle.HIGHEST_PROTOCOL)

    # Save an Xarray
    # with open(Path(dataPath, 'n2_3sg_0.1-50.1eV_A2_XR.pickle'), 'wb') as handle:
    dataFile = setTmpDataFile
    dataFile = dataFile +'_XR.pickle'

    print(dataFile)

    with open(dataFile, 'wb') as handle:
        pickle.dump(data, handle, protocol=pickle.HIGHEST_PROTOCOL)


def test_pickle_read(dataSingle, setTmpDataFile):

    data = dataSingle

    # Read back in to test
    # with open(Path(dataPath, 'n2_3sg_0.1-50.1eV_A2_XR.pickle'), 'rb') as handle:
    dataFile = setTmpDataFile
    dataFile = dataFile +'_XR.pickle'

    with open(dataFile, 'rb') as handle:
        dataDictPklIn = pickle.load(handle)

    assert data.equals(dataDictPklIn)  # This is True
    assert (dataDictPklIn - data).max() == 0   # And data looks OK.


def test_hdf5_write(dataSingle,setTmpDataFile):

    data = dataSingle

    dataFile = setTmpDataFile
    dataFile = dataFile +'_XR.h5'

    # ep.IO.writeXarray(data, fileName = 'n2_3sg_0.1-50.1eV_A2_XR.h5', filePath = None,
    #                     engine = 'hdf5')   # Default case set as: engine = 'h5netcdf', forceComplex = False
    ep.IO.writeXarray(data, fileName = dataFile,
                        engine = 'hdf5')   # Default case set as: engine = 'h5netcdf', forceComplex = False


def test_hdf5_read(dataSingle,setTmpDataFile):

    data = dataSingle

    dataFile = setTmpDataFile
    dataFile = dataFile +'_XR.h5'

    # HDF5 file reader, note this returns both dict and Xarray formats
    # dataInDict, dataInXR = ep.IO.readXarray(fileName = 'n2_3sg_0.1-50.1eV_A2_XR.h5',
    #                                         filePath = None, engine = 'hdf5')
    dataInDict, dataInXR = ep.IO.readXarray(fileName = dataFile,
                                            engine = 'hdf5')


    # assert data.equals(dataDictPklIn)  # This fails
    # assert data.dims == dataInXR.dims   # This may fail if dim ordering changed.
    assert not (set(dataInXR.dims) - set(data.dims))  # Use sets to ignore ordering, returns null if OK.
    assert (dataInXR - data).max() == 0   # And data looks OK.




def loadRef(fileName = None):
    """
    Print message to compare with readMatEle print output

    NOTE: not working, need to capture ref text more carefully.

    UPDATE: OK with file ref. (Write from calling function initially.)
    """
    if fileName is not None:
        with open(fileName, "r") as text_file:
        #     text_file.write(captured.out)
            ref = text_file.read()
    else:
        ref = "No file set"

    return ref

#     refPrint = """*** ePSproc readMatEle(): scanning files for DumpIdy segments.
#
# *** Scanning file(s)
# ['/home/femtolab/github/ePSproc/data/photoionization/n2_3sg_0.1-50.1eV_A2.inp.out']
#
# *** FileListSort
#   Prefix: /home/femtolab/github/ePSproc/data/photoionization/n2_3sg_0.1-50.1eV_A2.inp.out
#   1 groups.
#
# *** Reading ePS output file:  /home/femtolab/github/ePSproc/data/photoionization/n2_3sg_0.1-50.1eV_A2.inp.out
# *** IO.fileParse() found 1 segments with
# 	Start: ScatEng
# 	End: ['#'].
# Expecting 51 energy points.
# *** IO.fileParse() found 2 segments with
# 	Start: ScatSym
# 	End: ['FileName', '\n'].
# Expecting 2 symmetries.
# *** IO.fileParse() found 102 segments with
# 	Start: DumpIdy - dump
# 	End: ['+ Command', 'Time Now'].
# Found 102 dumpIdy segments (sets of matrix elements).
#
# Processing segments to Xarrays...
# Processed 102 sets of DumpIdy file segments, (0 blank)
# """

#     refPrint = r"""
# *** ePSproc readMatEle(): scanning files for DumpIdy segments.
#
# *** Scanning file(s)
# ['/home/femtolab/github/ePSproc/data/photoionization/n2_3sg_0.1-50.1eV_A2.inp.out']
#
# *** FileListSort
#   Prefix: /home/femtolab/github/ePSproc/data/photoionization/n2_3sg_0.1-50.1eV_A2.inp.out
#   1 groups.
#
# *** Reading ePS output file:  /home/femtolab/github/ePSproc/data/photoionization/n2_3sg_0.1-50.1eV_A2.inp.out
# *** IO.fileParse() found 1 segments with
# 	Start: ScatEng
# 	End: ['#'].
# Expecting 51 energy points.
# *** IO.fileParse() found 2 segments with
# 	Start: ScatSym
# 	End: ['FileName', '\n'].
# Expecting 2 symmetries.
# *** IO.fileParse() found 102 segments with
# 	Start: DumpIdy - dump
# 	End: ['+ Command', 'Time Now'].
# Found 102 dumpIdy segments (sets of matrix elements).
#
# Processing segments to Xarrays...
# Processed 102 sets of DumpIdy file segments, (0 blank)
# """

    # refPrint = """('*** ePSproc readMatEle(): scanning files for DumpIdy segments.\n'\n '\n'\n '*** Scanning file(s)\n'\n "['/home/femtolab/github/ePSproc/data/photoionization/n2_3sg_0.1-50.1eV_A2.inp.out']\n"\n '\n'\n '*** FileListSort \n'\n '  Prefix: '\n '/home/femtolab/github/ePSproc/data/photoionization/n2_3sg_0.1-50.1eV_A2.inp.out \n'\n '  1 groups.\n'\n '\n'\n '*** Reading ePS output file:  '\n '/home/femtolab/github/ePSproc/data/photoionization/n2_3sg_0.1-50.1eV_A2.inp.out\n'\n '*** IO.fileParse() found 1 segments with \n'\n '\tStart: ScatEng\n'\n "\tEnd: ['#'].\n"\n 'Expecting 51 energy points.\n'\n '*** IO.fileParse() found 2 segments with \n'\n '\tStart: ScatSym\n'\n "\tEnd: ['FileName', '\\n'].\n"\n 'Expecting 2 symmetries.\n'\n '*** IO.fileParse() found 102 segments with \n'\n '\tStart: DumpIdy - dump\n'\n "\tEnd: ['+ Command', 'Time Now'].\n"\n 'Found 102 dumpIdy segments (sets of matrix elements).\n'\n '\n'\n 'Processing segments to Xarrays...\n'\n 'Processed 102 sets of DumpIdy file segments, (0 blank)\n')"""
    #
    # return refPrint
