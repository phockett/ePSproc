"""
ePSproc utility functions.

Set of tools for assignment, sorting, normalisation and conversion.

16/03/20    Converted to submodule, mainly split out from old util.py, plus some new functions.
            Imports may be buggy...
14/10/19    Added string replacement function (generic)
11/08/19    Added matEleSelector

"""

# import numpy as np
# import re
# import scipy.constants

# # Package fns.
# from epsproc.basicPlotters import molPlot

# Set explict imports
# from .conversion import multiDimXrToPD, conv_ev_atm, conv_ev_nm
# from .listFuncs import matEdimList, BLMdimList, eulerDimList, ADMdimList, dataTypesList, genLM
# from .summary import jobSummary, lmSymSummary
# from .selectors import matEleSelector, dataGroupSel
# from .misc import stringRepMap, arraySort2D

# Autogen with mkinit, https://github.com/Erotemic/mkinit
# <AUTOGEN_INIT>
from epsproc.util import conversion
from epsproc.util import listFuncs
from epsproc.util import misc
from epsproc.util import selectors
from epsproc.util import summary

from epsproc.util.conversion import (conv_ev_atm, conv_ev_nm, multiDimXrToPD,orb3DCoordConv)
from epsproc.util.listFuncs import (ADMdimList, BLMdimList, dataTypesList,
                                    eulerDimList, genLM, matEdimList,YLMtype)
from epsproc.util.misc import (arraySort2D, stringRepMap,)
from epsproc.util.selectors import (dataGroupSel, matEleSelector,)
from epsproc.util.summary import (jobSummary, lmSymSummary,)

__all__ = ['ADMdimList', 'BLMdimList', 'arraySort2D', 'conv_ev_atm',
           'conv_ev_nm', 'conversion', 'dataGroupSel', 'dataTypesList',
           'eulerDimList', 'genLM', 'jobSummary', 'listFuncs', 'lmSymSummary',
           'matEdimList', 'matEleSelector', 'misc', 'multiDimXrToPD',
           'orb3DCoordConv',
           'selectors', 'stringRepMap', 'summary']
# </AUTOGEN_INIT>
