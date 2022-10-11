# import os
# dataPath = os.path.join(sys.path[-1], 'data', 'photoionization', 'n2_multiorb')

import pytest
from pathlib import Path

# ePSproc
import epsproc as ep
from epsproc.classes.base import ePSbase

@pytest.fixture
def setDataPath():
    # Set data path
    # Note this is set here from ep.__path__, but may not be correct in all cases - depends on where the Github repo is.
    epDemoDataPath = Path(ep.__path__[0]).parent/'data'

    # dataPath = os.path.join(epDemoDataPath, 'photoionization', 'n2_multiorb')
    dataPath = Path(epDemoDataPath, 'photoionization', 'n2_multiorb')

    return dataPath


def test_base_init(setDataPath):

    data = ePSbase(setDataPath, verbose = 1)

    data.scanFiles()
