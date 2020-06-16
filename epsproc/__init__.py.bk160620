__version__ = '1.2.5-dev'

# Import shared packages - actually, should be at module (file) level?
# https://stackoverflow.com/questions/8165703/python-imports-across-modules-and-global-variables
# OR GLOBAL...?
# https://stackoverflow.com/questions/11990556/how-to-make-global-imports-from-a-function
#def my_imports(module_name):
#    globals()[module_name] = __import__(module_name)

# import numpy as np


# Import functions
#from epsproc import ePSproc_IO as epIO
#from epsproc import ePSproc_util as epUT
#from epsproc import ePSproc_sphCalc as epSC
#from epsproc import ePSproc_sphPlot as epSP
#from epsproc import ePSproc_MFPAD as epMF

#import epsproc.ePSproc_IO as IO
#import epsproc.ePSproc_util as UT
#import epsproc.ePSproc_sphCalc as SC
#import epsproc.ePSproc_sphPlot as SP
#import epsproc.ePSproc_MFPAD as MF

# Direct * import for modules - probably a bad idea.
# from epsproc.IO import *
# from epsproc.util import *
# from epsproc.basicPlotters import *
# from epsproc.sphCalc import *
# from epsproc.sphPlot import *
# from epsproc.MFPAD import *
# from epsproc.MFBLM import *
# from epsproc.AFBLM import afblm, AFBLMCalcLoop
# from epsproc.conversion import multiDimXrToPD



# TODO - move from individual modules.
# from epsproc.util import matEleSelector

# Testing mkinit, https://github.com/Erotemic/mkinit
# On Stimpy, run with
# > mkinit D:\code\github\ePSproc\epsproc > mkinitOut.txt
# Writing directly to __init__.py sets autogen code only?
#
# <AUTOGEN_INIT>
from epsproc import AFBLM
from epsproc import IO
from epsproc import MFBLM
from epsproc import MFPAD
from epsproc import basicPlotters
from epsproc import geomFunc
from epsproc import sphCalc
from epsproc import sphPlot
from epsproc import util

from epsproc.AFBLM import (AFBLMCalcLoop, Wigner3jCached,
                           Wigner_D_element_Cached, afblm, blmXarray,)
from epsproc.IO import (EDCSFileParse, EDCSSegParse, EDCSSegsParseX,
                        dumpIdyFileParse, dumpIdySegParse, dumpIdySegsParseX,
                        fileParse, getCroFileParse, getCroSegParse,
                        getCroSegsParseX, getFiles, headerFileParse,
                        matEleGroupDim, matEleGroupDimX, matEleGroupDimXnested,
                        molInfoParse, parseLineDigits, readMatEle, readOrb3D,
                        readOrbCoords, readOrbData, readOrbElements,
                        readOrbHeader, readXarray, scatEngFileParse,
                        symFileParse, writeOrb3Dvtk, writeXarray,)
from epsproc.MFBLM import (MFBLMCalcLoop, Wigner3jCached,
                           Wigner_D_element_Cached, blmXarray, mfblm,
                           mfblmEuler,)
from epsproc.MFPAD import (mfpad,)
from epsproc.basicPlotters import (Arrow3D, BLMplot, lmPlot, molPlot,
                                   symListGen,)
from epsproc.sphCalc import (TKQarrayRot, TKQarrayRotX, setADMs, setPolGeoms,
                             sphCalc, wDcalc,)
from epsproc.sphPlot import (plotTypeSelector, sphFromBLMPlot, sphPlotHV,
                             sphPlotMPL, sphPlotPL, sphSumPlotX, sphToCart,)

__all__ = ['AFBLM', 'AFBLMCalcLoop', 'Arrow3D', 'BLMplot', 'EDCSFileParse',
           'EDCSSegParse', 'EDCSSegsParseX', 'IO', 'MFBLM', 'MFBLMCalcLoop',
           'MFPAD', 'TKQarrayRot', 'TKQarrayRotX', 'Wigner3jCached',
           'Wigner3jCached', 'Wigner_D_element_Cached',
           'Wigner_D_element_Cached', 'afblm', 'basicPlotters', 'blmXarray',
           'blmXarray', 'dumpIdyFileParse', 'dumpIdySegParse',
           'dumpIdySegsParseX', 'fileParse', 'geomFunc', 'getCroFileParse',
           'getCroSegParse', 'getCroSegsParseX', 'getFiles', 'headerFileParse',
           'lmPlot', 'matEleGroupDim', 'matEleGroupDimX',
           'matEleGroupDimXnested', 'mfblm', 'mfblmEuler', 'mfpad',
           'molInfoParse', 'molPlot', 'parseLineDigits', 'plotTypeSelector',
           'readMatEle', 'readOrb3D', 'readOrbCoords', 'readOrbData',
           'readOrbElements', 'readOrbHeader', 'readXarray',
           'scatEngFileParse', 'setADMs', 'setPolGeoms', 'sphCalc', 'sphCalc',
           'sphFromBLMPlot', 'sphPlot', 'sphPlotHV', 'sphPlotMPL', 'sphPlotPL',
           'sphSumPlotX', 'sphToCart', 'symFileParse', 'symListGen', 'util',
           'wDcalc', 'writeOrb3Dvtk', 'writeXarray']
# </AUTOGEN_INIT>
