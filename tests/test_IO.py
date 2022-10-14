# Basic IO test routines.
# Mainly from https://epsproc.readthedocs.io/en/dev/dataStructures/ePSproc_dataStructures_IO_demo_280622.html

import pytest
from pathlib import Path

# ePSproc
import epsproc as ep

# NOTE - should centralise this routine!
# @pytest.fixture(scope="session")
@pytest.fixture(scope="module")  # Scoping here doesn't help with class issues below.
def setDataPath():
    # Set data path
    # Note this is set here from ep.__path__, but may not be correct in all cases - depends on where the Github repo is.
    epDemoDataPath = Path(ep.__path__[0]).parent/'data'

    # dataPath = os.path.join(epDemoDataPath, 'photoionization', 'n2_multiorb')
    dataPath = Path(epDemoDataPath, 'photoionization')

    return dataPath


def test_ePS_read(setDataPath, capsys):
    # Load data from modPath\data
    dataPath = setDataPath
    dataFile = Path(dataPath, 'n2_3sg_0.1-50.1eV_A2.inp.out')  # Set for sample N2 data for testing

    # Scan data file
    dataSet = ep.readMatEle(fileIn = dataFile.as_posix())
    data = dataSet[0]

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

    assert captured.out == loadRef(fileName = 'ePS_read_ref.txt')
    # print('OK')


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
