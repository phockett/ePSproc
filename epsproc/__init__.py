__version__ = '1.1.0'

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
from epsproc.ePSproc_IO import *
from epsproc.ePSproc_util import *
from epsproc.ePSproc_sphCalc import *
from epsproc.ePSproc_sphPlot import *
from epsproc.ePSproc_MFPAD import *



# TODO - move from individual modules.
# from epsproc.ePSproc_util import matEleSelector