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
from epsproc.geomFunc import geomCalc
from epsproc.geomFunc import geomUtils
from epsproc.geomFunc import mfblmGeom
# from epsproc.geomFunc import mfblmGeom_Dev
from epsproc.geomFunc import w3jVecMethods

from epsproc.geomFunc.geomCalc import (EPR, MFproj, betaTerm, remapllpL,
                                       setPhaseConventions, w3jTable,)
from epsproc.geomFunc.geomUtils import (genllL, genllpMatE, selQNsRow,)
from epsproc.geomFunc.mfblmGeom import (mfblmXprod,)
from epsproc.geomFunc.mfblmGeom_Dev import (mfblmXprod,)
from epsproc.geomFunc.w3jVecMethods import (Wigner3jQNs, w3jguVecCPU,
                                            w3jprange,)

__all__ = ['EPR', 'MFproj', 'Wigner3jQNs', 'betaTerm', 'genllL', 'genllpMatE',
           'geomCalc', 'geomUtils', 'mfblmGeom', 'mfblmGeom_Dev', 'mfblmXprod',
           'mfblmXprod', 'remapllpL', 'selQNsRow', 'setPhaseConventions',
           'w3jTable', 'w3jVecMethods', 'w3jguVecCPU', 'w3jprange']
# </AUTOGEN_INIT>
