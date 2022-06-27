
__version__ = '1.3.2-dev'


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

# Autogen with mkinit, https://github.com/Erotemic/mkinit
# <AUTOGEN_INIT>

from epsproc import AFBLM
from epsproc import IO
from epsproc import MFBLM
from epsproc import MFPAD
from epsproc import basicPlotters
from epsproc import efield
from epsproc import geomFunc
from epsproc import sphCalc
from epsproc import sphPlot
from epsproc import util
from epsproc import vol

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
from epsproc.efield import (Efield, efields,)
from epsproc.geomFunc import (EPR, MFproj, Wigner3jQNs, afblmGeom, afblmXprod,
                              betaTerm, deltaLMKQS, genKQSterms,
                              genKQStermsFromTensors, genllL, genllLList,
                              genllpMatE, geomCalc, geomUtils, mfblmGeom,
                              mfblmXprod, mfblmXprod, remapllpL,
                              selQNsRow, setPhaseConventions,
                              setPhaseConventions, w3jTable, w3jVecMethods,
                              w3jguVecCPU, w3jprange,)
from epsproc.sphCalc import (TKQarrayRot, TKQarrayRotX, setADMs, setPolGeoms,
                             sphCalc, wDcalc,)
from epsproc.sphPlot import (plotTypeSelector, sphFromBLMPlot, sphPlotHV,
                             sphPlotMPL, sphPlotPL, sphSumPlotX, sphToCart,)
from epsproc.util import (ADMdimList, BLMdimList, arraySort2D, conv_ev_atm,
                          conv_ev_nm, conversion, dataGroupSel, dataTypesList,
                          eulerDimList, genLM, jobSummary, listFuncs,
                          lmSymSummary, matEdimList, matEleSelector, misc,
                          multiDimXrToPD, selectors, stringRepMap, summary,)

__all__ = ['ADMdimList', 'AFBLM', 'AFBLMCalcLoop', 'Arrow3D', 'BLMdimList',
           'BLMplot', 'EDCSFileParse', 'EDCSSegParse', 'EDCSSegsParseX', 'EPR',
           'Efield', 'IO', 'MFBLM', 'MFBLMCalcLoop', 'MFPAD', 'MFproj',
           'TKQarrayRot', 'TKQarrayRotX', 'Wigner3jCached', 'Wigner3jCached',
           'Wigner3jQNs', 'Wigner_D_element_Cached', 'Wigner_D_element_Cached',
           'afblm', 'afblmGeom', 'afblmXprod', 'arraySort2D', 'basicPlotters',
           'betaTerm', 'blmXarray', 'blmXarray', 'conv_ev_atm', 'conv_ev_nm',
           'conversion', 'dataGroupSel', 'dataTypesList', 'deltaLMKQS',
           'dumpIdyFileParse', 'dumpIdySegParse', 'dumpIdySegsParseX',
           'efield', 'efields', 'eulerDimList', 'fileParse', 'genKQSterms',
           'genKQStermsFromTensors', 'genLM', 'genllL', 'genllLList',
           'genllpMatE', 'geomCalc', 'geomFunc', 'geomUtils',
           'getCroFileParse', 'getCroSegParse', 'getCroSegsParseX', 'getFiles',
           'headerFileParse', 'jobSummary', 'listFuncs', 'lmPlot',
           'lmSymSummary', 'matEdimList', 'matEleGroupDim', 'matEleGroupDimX',
           'matEleGroupDimXnested', 'matEleSelector', 'mfblm', 'mfblmEuler',
           'mfblmGeom', 'mfblmGeom_Dev', 'mfblmXprod', 'mfblmXprod', 'mfpad',
           'misc', 'molInfoParse', 'molPlot', 'multiDimXrToPD',
           'parseLineDigits', 'plotTypeSelector', 'readMatEle', 'readOrb3D',
           'readOrbCoords', 'readOrbData', 'readOrbElements', 'readOrbHeader',
           'readXarray', 'remapllpL', 'scatEngFileParse', 'selQNsRow',
           'selectors', 'setADMs', 'setPhaseConventions',
           'setPhaseConventions', 'setPolGeoms', 'sphCalc', 'sphCalc',
           'sphFromBLMPlot', 'sphPlot', 'sphPlotHV', 'sphPlotMPL', 'sphPlotPL',
           'sphSumPlotX', 'sphToCart', 'stringRepMap', 'summary',
           'symFileParse', 'symListGen', 'util', 'vol', 'w3jTable',
           'w3jVecMethods', 'w3jguVecCPU', 'w3jprange', 'wDcalc',
           'writeOrb3Dvtk', 'writeXarray']

# </AUTOGEN_INIT>
