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
    df = pd.read_csv(fileIn, skiprows=lineStart[0]-1, # skiprows seems more consistent than 'header'
                     sep='\t+',   # use reg ex to allow for inconsistent tabbing
                     engine='python',  # Set engine to allow reg ex without warning
                     index_col=cols)  # Add index
    
    # Read file header
    with open(fileIn) as fp:
        # lines = fp.readlines(100)
        head = [next(fp) for _ in range(lineStop[0])]
        
        # Get total number of (data)lines post header
        lines = len(fp.readlines())
        
    # TODO - also check file length?
        
    # Set metadata
    df.attrs['header']=head
    df.attrs['file']=fileIn
    df.attrs['headerlines']=[lineStart, lineStop]
    df.attrs['lines']=lines
    
    return df