"""
ePSproc Gamma functions: IO

- Read legacy Gamma files.

15/07/24

"""

# Imports
import os
import pandas as pd

from epsproc import fileParse


def readGamma(fileName, filePath = None, 
              headerLines=None, cols=None):
    """
    Read gamma-parameters from file, and convert to Pandas DataFrame.
    
    Legacy file format:
    
        nh3 ion simulation gamma values
        Fri Jul 10 16:35:47 2009

        l=0:2:6
        |N,K>: X(0,0) -> B(2,1) -> X+(1,0)

        gamma			L	M	l1	lambda1	ml1	l2	lambda2	ml2	betaTerm

        0.074074		0	0	0	0	0	0	0	0	1.000000
        0.015923		2	0	0	0	0	2	-2	0	1.000000
        ....
    
    
    Parameters
    ----------
    fileName : str or path object
        Input filename.
        
    filePath : str or path object, optional, default = None
        Full path.
        If not set read from current dir.
        
    headLines : int, optional, default = None
        Header lines to skip at start.
        If None will be tested and set automatically.
        
    cols : list, optional, default = None
        Columns from file to use as DataFrame index
        Default case uses: cols = list(range(1,9))
       
       
    Returns
    -------
    DataFrame
        Results in dataframe.
        See df.attrs for additional file info
        
    
    TODO:
    - Path() for files
    - Logging/info to add
    
    """
    
    
    if filePath is not None:
        # filePath = os.getcwd()
        fileName = os.path.join(filePath, fileName)

    if not os.path.exists(fileName):
        print(f"*** No file found for {fileName}, skipping file read.")
        return None

    # Check header lines using existing routine.
    # TODO: should just use RE parsing here...?
    if headerLines is None:
        # ([lineStart, lineStop], segments) = checkHeader
        ([lineStart, lineStop], segments) = fileParse(fileName,
                                                      startPhrase = 'gamma',
                                                      endPhrase = '\n')
        
    if cols is None:
        cols = list(range(1,9))
    
    # Read CSV with PD
    df = pd.read_csv(fileName, skiprows=lineStart[0]-1, # skiprows seems more consistent than 'header'
                     sep='\t+',   # use reg ex to allow for inconsistent tabbing
                     engine='python',  # Set engine to allow reg ex without warning
                     index_col=cols)  # Add index
    
    # Read file header
    with open(fileName) as fp:
        # lines = fp.readlines(100)
        head = [next(fp) for _ in range(lineStop[0])]
        
        # Get total number of (data)lines post header
        lines = len(fp.readlines())
        
    # TODO - also check file length?
        
    # Set metadata
    df.attrs['header']=head
    df.attrs['file']=fileName
    df.attrs['headerlines']=[lineStart, lineStop]
    df.attrs['lines']=lines
    df.attrs['dataType']="gamma"
    df.attrs['source']="file"    # 09/08/24 Set this to allow quick switch on legacy gamma from file vs. new python calcs in ancillary functions.
    df.attrs['legacyGamma']=True
    
    return df


def toePSprocClass(dataIn, dimMap={'L':'l','M':'m','variable':'t'}, 
                   conformDims = True, dropDims = None,
                   key = None, dataType = None, **kwargs):
    """
    Convert gamma and derivatives from stand-alone PD dataframe to ePSproc class data object (with data in Xarray).
    
    TODO: see PEMtk functionality for additional methods.
    
    Parameters
    ----------
    dataIn : PD DataFrame
        Gamma data, or derivatives, as Pandas DataFrame
        
    dimMap : dict, optional, default = {'L':'l','M':'m','variable':'t'}
        Dim remapping for conversion.
        Default for AFBLM case.
        
    conformDims : bool, optional, default = True
        Force XR dims to match dataType specification if True, via call to :py:func:`pemtk.sym._util.toePSproc`.
        
    dropDims : str or list, optional, default = None
        If passed, run data.drop_vars(dropDims)
        
    key : str, optional, default = None
        Key for output data in ePSproc class.
        If None, use default key='gamma'
        
    dataType : str, optional, default = None
        dataType to set for output.
        If None try and use dataXR.attrs['dataType'], or default to 'AFBLM' if not set.
        
    Returns
    -------
    dataOut : epsproc.classes.multiJob
        Class with data set as `dataOut.data[key][dataType]`.
    
    """
    
    # Convert to XR
    dataXR = dataIn.to_xarray().to_array()
    
    # Set dataType for class if not preset.
    # NOTE this may be missing in some cases, set AFBLM as default
    if dataType is None:
        try:
            dataType = dataXR.attrs['dataType']
        except KeyError:
            print(f"*** Warning: dataIn.attrs['dataType'] not found, setting dataType='AFBLM', or pass dataType to override.")
            dataType = 'AFBLM'
    
    # Test conversion with existing functionality...
    from epsproc.util.conversion import multiDimXrFromDict
    # multiDimXrFromDict(betaOutNorm.to_dict())  # Needs XR dict format

    if conformDims:
        # This works with some effort
        from pemtk.sym._util import toePSproc
        # coeffs = {'XR':betaXR}
        dataXRremapped = toePSproc({'XR':dataXR}, dimMap=dimMap, dataType=dataType, **kwargs)
    
    else:
        dataXRremapped = dataXR
        
    if dropDims is not None:
        dataXRremapped = dataXRremapped.drop_vars(dropDims)
    
    # Push to ePSproc class
    from epsproc.classes.multiJob import ePSmultiJob
    dataOut = ePSmultiJob()
    
    if key is None:
        key = 'gamma'
    
    dataOut.data[key] = {dataType: dataXRremapped}  #.drop_vars('Euler')}
    
    return dataOut

    # data.BLMplot(xDim='t', backend='hv', hvType='line') #, addADMs=False)  #, col='Eke')