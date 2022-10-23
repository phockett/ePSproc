import pytest
from pathlib import Path
import epsproc as ep

# @pytest.fixture(scope="module")  # Scoping here doesn't help with class issues below.
@pytest.fixture(scope="session")
def setDataPath():
    # Set data path
    # Note this is set here from ep.__path__, but may not be correct in all cases - depends on where the Github repo is.
    epDemoDataPath = Path(ep.__path__[0]).parent/'data'

    # dataPath = os.path.join(epDemoDataPath, 'photoionization', 'n2_multiorb')
    dataPath = Path(epDemoDataPath, 'photoionization')

    return dataPath

# @pytest.fixture(scope="module")
@pytest.fixture(scope="session")
def setDataFileName():
    return 'n2_3sg_0.1-50.1eV_A2'

# @pytest.fixture(scope="module")
@pytest.fixture(scope="session")
def setDataFile(setDataPath,setDataFileName):
    dataPath = setDataPath
    return Path(dataPath,setDataFileName).as_posix()  # Return POSIX to handle multiple suffixes later.
                                                             # Or use .with_name?

# Use tmp_path fixture for test data files, see https://docs.pytest.org/en/7.1.x/how-to/tmp_path.html
# tmp_path_factory should be better for persistence? https://docs.pytest.org/en/7.1.x/how-to/tmp_path.html#the-tmp-path-factory-fixture

# @pytest.fixture()
# def setTmpDataFile(tmp_path,setDataFileName,tmp_path_factory):
@pytest.fixture(scope="session")
def setTmpDataFile(setDataFileName,tmp_path_factory):
    # dataPath = tmp_path  # THIS WORKS, but note that this also sets tmp dirs PER TEST FUNCTION.
    dataPath = tmp_path_factory.mktemp("data") / setDataFileName  # SAME ISSUE - sets DATAN per N run.
                                                                  # NEEDS scope = "session" to use single dir over all tests.
                                                                  # Note this may fail if other fixtures are NOT session scoped.
    return dataPath.as_posix()

    # return Path(dataPath,setDataFileName).as_posix()  # Return POSIX to handle multiple suffixes later.
                                                             # Or use .with_name?

# @pytest.fixture(scope="session")
# def image_file(tmp_path_factory):
#     img = compute_expensive_image()
#     fn = tmp_path_factory.mktemp("data") / "img.png"
#     img.save(fn)
#     return fn


# Data as fixutre?
# @pytest.fixture(scope="module")
@pytest.fixture(scope="session")
# def dataSingle(setDataPath):
def dataSingle(setDataFile):
    # dataPath = setDataPath
    # dataFile = Path(dataPath, 'n2_3sg_0.1-50.1eV_A2.inp.out')  # Set for sample N2 data for testing

    # dataFile = setDataFile.with_suffix('.inp.out')
    dataFile = setDataFile+'.inp.out'

    # Scan data file
    # dataSet = ep.readMatEle(fileIn = dataFile.as_posix())
    dataSet = ep.readMatEle(fileIn = dataFile)
    data = dataSet[0]

    return data
