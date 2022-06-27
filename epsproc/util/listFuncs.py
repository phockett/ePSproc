"""
ePSproc list functions.

Define data types and generate lists.

Main function for data types is :py:func:`dataTypesList`

"""

import numpy as np
# import ..geomFunc.setPhaseConventions
from epsproc.geomFunc import setPhaseConventions

#************************** LIST functions

# Return list of standard dataArray dims for matrix elements
def matEdimList(sType = 'stacked'):
    """
    Return standard list of dimensions for matrix elements.

    Parameters
    ----------
    sType : string, optional, default = 'stacked'
        Selected 'stacked' or 'unstacked' dimensions.
        Set 'sDict' to return a dictionary of unstacked <> stacked dims mappings for use with xr.stack({dim mapping}).

    Returns
    -------
    list : set of dimension labels.

    """
    if sType == 'stacked':
        # stackedDims
        return ['LM', 'Eke', 'Sym', 'mu', 'it', 'Type']

    elif sType == 'sDict':
        return {'LM':['l','m'],'Sym':['Cont', 'Targ', 'Total']}

    else:
        # unstackedDims
        return ['l','m', 'Eke', 'Cont', 'Targ', 'Total', 'mu', 'it', 'Type']

# Return list of standard dataArray dims for BLM values
def BLMdimList(sType = 'stacked'):
    """
    Return standard list of dimensions for calculated BLM.

    Parameters
    ----------
    sType : string, optional, default = 'stacked'
        Selected 'stacked' or 'unstacked' dimensions.
        Set 'sDict' to return a dictionary of unstacked <> stacked dims mappings for use with xr.stack({dim mapping}).

    Returns
    -------
    list : set of dimension labels.

    """
    if sType == 'stacked':
        # stackedDims
        return ['Euler', 'Eke', 'BLM']

    elif sType == 'sDict':
        return {'BLM':['l','m'],'Euler':['P','T','C']}

    else:
        # unstackedDims
        return ['Eke', 'l', 'm', 'P', 'T', 'C']

# Return list of standard dataArray dims for Euer angles
def eulerDimList(sType = 'stacked'):
    """
    Return standard list of dimensions for frame definitions, from :py:func:`epsproc.sphCalc.setPolGeoms()`.

    Parameters
    ----------
    sType : string, optional, default = 'stacked'
        Selected 'stacked' or 'unstacked' dimensions.
        Set 'sDict' to return a dictionary of unstacked <> stacked dims mappings for use with xr.stack({dim mapping}).

    Returns
    -------
    list : set of dimension labels.

    """
    if sType == 'stacked':
        # stackedDims
        return ['Euler']

    elif sType == 'sDict':
        return {'Euler':eulerDimList(sType = 'unstacked')}

    else:
        # unstackedDims
        return ['Label','P','T','C']

# Return list of standard dataArray dims for Euer angles
def ADMdimList(sType = 'stacked'):
    """
    Return standard list of dimensions for frame definitions, from :py:func:`epsproc.sphCalc.setADMs()`.

    Parameters
    ----------
    sType : string, optional, default = 'stacked'
        Selected 'stacked' or 'unstacked' dimensions.
        Set 'sDict' to return a dictionary of unstacked <> stacked dims mappings for use with xr.stack({dim mapping}).

    Returns
    -------
    list : set of dimension labels.

    """
    if sType == 'stacked':
        # stackedDims
        return ['ADM','t']

    elif sType == 'sDict':
        return {'ADM':['K','Q','S'],'t':'t'}

    else:
        # unstackedDims
        return ['K','Q','S','t']



#**************** Master list fn.
# Return a dict of allowed dataTypes, corresponding to ePS commands/file segments, or epsproc processed data.
def dataTypesList():
    """
    Return a dict of allowed dataTypes, corresponding to epsproc processed data.

    Each dataType lists 'source', 'desc' and 'recordType' fields.

    - 'source' fields correspond to ePS functions which get or generate the data.
    - 'desc' brief description of the dataType.
    - 'recordType' gives the required segment in ePS files (and associated parser). If the segment is not present in the source file, then the dataType will not be available.
    - 'def' provides definition function handle (if applicable).
    - 'dims' lists results from def(sType = 'sDict')

    TODO: best choice of data structure here?  Currently nested dictionary.

    """

    dataDict = {'BLM' :
                    {'source':'epsproc.MFBLM',
                    'desc':'Calcualted MF beta parameters from epsproc.MFBLM(), based on dipole matrix elements from ePS.',
                    'recordType':'DumpIdy',
                    'dims': BLMdimList(sType = 'sDict'),
                    'def': BLMdimList
                    },
                'matE' :
                    {'source':'epsproc.IO.readMatEle',
                    'desc':'Raw photoionization matrix elements from ePS, DumpIdy command and file segments.',
                    'recordType':'DumpIdy',
                    'dims': matEdimList(sType = 'sDict'),
                    'def': matEdimList
                    },
                'EDCS' :
                    {'source':'epsproc.IO.readMatEle',
                    'desc':'Electron scattering DCS results from ePS, EDCS command and file segments.',
                    'recordType':'EDCS'
                    },
                'XSect' :
                    {'source':'CrossSection',
                    'desc':'Photoionziation cross section results from ePS, GetCro command and CrossSection file segments.',
                    'recordType':'CrossSection'
                    },
                'TX' :
                    {'source':'epsproc.MFPAD.mfpad',
                    'desc':'MF continuum wavefunction (E,theta,phi), expanded or summed over (l,m). Abs(TX)^2 gives the MFPAD.',
                    'recordType':'DumpIdy'
                    },
                'Wigner' :
                    {'source':'epsproc.MFPAD.mfWignerDelay',
                    'desc':'Phase derivative of the MF continuum wavefunction (E,theta,phi), expanded or summed over (l,m), converted to temporal units.',
                    'recordType':'DumpIdy'
                    },
                'Euler' :
                    {'source':'epsproc.sphCalc.setPolGeoms()',
                     'desc':'Frame rotation definitions in terms of Euler angles and corresponding Quaternions.',
                     'dims':eulerDimList(sType = 'sDict'),
                     'def': eulerDimList
                     },
                 'ADM' :
                     {'source':'epsproc.sphCalc.setADMs()',
                      'desc':'Defines ADMs, A(K,Q,S), for aligned frame distributions.',
                      'dims':ADMdimList(sType = 'sDict'),
                      'def':ADMdimList
                      },
                'phaseCons' :
                     {'source':'epsproc.geomFunc.setPhaseConventions()',
                      'desc':'Defines sets of phase choices for calculations, implemented on in geomFuncs.',
                      'dims': setPhaseConventions(typeList = True),   # Currently supported types, should integrate in fn.
                      'defns': [setPhaseConventions(item) for item in setPhaseConventions(typeList = True)],
                      'def': setPhaseConventions
                      }
                }

    return dataDict

# Get ref dims for self
def getRefDims(data = None, refType = None, sType = 'sDict'):
    """
    Get ref dims from dataTypesList() using data or refType.

    Convenience wrapper for dataTypesList()[dataType]['def'](sType=sType).

    Parameters
    ----------

    data : optional, Xarray
        Data to use for dataType = data.attrs['dataType']

    refDims : optional, str, list or dict. Default = None.
        Reference dimensions
        - If None, use data.attrs['dataType']
        - If string, use dataTypesList()[refDims]['def'](sType='sDict')
        - Other types will be ignored and returned.

    """

    if refType is None:
        try:
            return dataTypesList()[data.attrs['dataType']]['def'](sType=sType)
        except:
            print("*** Input Xarray missing self.attrs['dataType'], can't get reference dimensions. To fix set missing attr, or pass refDims as string, list or dict.")
            return 0

    if isinstance(refType, str):
        return dataTypesList()[refType]['def'](sType=sType)

    # For other refTypes just return.
    # This allows for ignoring dictionary or list refDim cases.
    else:
        return refType


#****** QN lists

# Generate LM index from 0:Lmax
# Should be able to use sf.LM_range here, but throwing type errors for some reason.
def genLM(Lmax):
    """
    Return array of (L,M) up to supplied Lmax

    TODO: add return type options, include conversion to SHtools types.
    """

    # Loop version
#    LM = np.empty((Lmax * (Lmax + 2) + 1, 2), dtype=int)
#    n=0
#    for L in np.arange(Lmax+1):
#        for M in np.arange(-L,L+1):
#            LM[n,:] = [L,M]
#            n+=1

    # List version
    LM = []
    for L in np.arange(Lmax+1):
        LM.extend(np.c_[np.tile(L,2*L+1), np.arange(-L,L+1)])
#        LM.append([np.tile(L,2*L+1), np.arange(-L,L+1)])

    return np.asarray(LM)
