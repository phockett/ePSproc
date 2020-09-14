"""
ePSproc geomFunc

Codes for geometric functions:

- Wigner 3j
- Wigner D
- Spherical harmonics
- EPR tensors
- BLM tensors
etc.

Status
------
In progress Feb. 2020 for MFBLM development. Will eventually supercede/subsume existing functions in sphCalc.py.


"""


# Tags for mkinit, https://github.com/Erotemic/mkinit
# <AUTOGEN_INIT>
from epsproc.geomFunc import afblmGeom
from epsproc.geomFunc import geomCalc
from epsproc.geomFunc import geomUtils
from epsproc.geomFunc import mfblmGeom
from epsproc.geomFunc import mfblmGeom_Dev
# from epsproc.geomFunc import phaseCons-verified_130420
from epsproc.geomFunc import w3jVecMethods

from epsproc.geomFunc.afblmGeom import (afblmXprod,)
from epsproc.geomFunc.geomCalc import (EPR, MFproj, betaTerm, deltaLMKQS,
                                       remapllpL, setPhaseConventions,
                                       w3jTable,)
from epsproc.geomFunc.geomUtils import (genKQSterms, genKQStermsFromTensors,
                                        genllL, genllLList, genllpMatE,
                                        selQNsRow,)
from epsproc.geomFunc.mfblmGeom import (mfblmXprod,)
from epsproc.geomFunc.mfblmGeom_Dev import (mfblmXprod,)
# from epsproc.geomFunc.phaseCons-verified_130420 import (setPhaseConventions,)
from epsproc.geomFunc.w3jVecMethods import (Wigner3jQNs, w3jguVecCPU,
                                            w3jprange,)

__all__ = ['EPR', 'MFproj', 'Wigner3jQNs', 'afblmGeom', 'afblmXprod',
           'betaTerm', 'deltaLMKQS', 'genKQSterms', 'genKQStermsFromTensors',
           'genllL', 'genllLList', 'genllpMatE', 'geomCalc', 'geomUtils',
           'mfblmGeom', 'mfblmGeom_Dev', 'mfblmXprod', 'mfblmXprod',
           'remapllpL', 'selQNsRow',
           'setPhaseConventions', 'setPhaseConventions', 'w3jTable',
           'w3jVecMethods', 'w3jguVecCPU', 'w3jprange']
# </AUTOGEN_INIT>
