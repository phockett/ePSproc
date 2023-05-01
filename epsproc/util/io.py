"""
ePSproc IO util functions

Various tools for use in file IO.

01/05/23    Added getFilesFromURLs.
            See also epsdata module for Zenodo-specific downloader.

27/06/22    Split out from core IO.py, to extend backend options and support.
            Now additionally wrapped therein for flexible handling of multiple backends.

"""

from datetime import datetime as dt # Import datetime.datetime for now() function
import wget
from pathlib import Path

from epsproc.util.epsdata import convert_bytes


def setTimeStampedFileName(outStem = None, n = None, ext = None, timeString = None, timeFormat = '%Y-%m-%d_%H-%M-%S'):
    """
    Set unique filename as f'{outStem}_n{n}_{timeString.strftime("%d%m%y_%H-%M-%S")}.{ext}'

    Parameters
    ----------

    outStem : str, optional, default = None
        Stem for output file.
        If None, set to 'ep'

    n : int, optional, default = None
        Int index to include in file name.
        If None, this will be omitted.

    ext : str, optional, default = None
        File ending.
        If None, this will be ommitted

    timeString : Datatime object, default = None
        Timestamp for the file.
        If None, current time will be used.

    timeFormat : Datatime format string, optional, default = '%Y-%m-%d_%H-%M-%S'

    """

    if outStem is None:
        outStem = 'ep'

    if timeString is None:
        timeString = dt.now()

    # Use n for multiple file checkpointing - may want a more sophisticated routine for this however.
    if n is not None:
        fName = f'{outStem}_n{n}_{timeString.strftime(timeFormat)}'
    else:
        fName = f'{outStem}_{timeString.strftime(timeFormat)}'

    if ext is not None:
        fName = fName + f'.{ext}'

    return fName



def getFilesFromURLs(urls, dataName=None, dataPath=None, ghRaw = True):
    """
    Pull (dataset) files from dictionary of URLs.

    Note this adds '?raw=true' if ghRaw=True, and 'github.com' is in the URL.

    Parameters
    ----------
    urls : list or dict
        URLs to pull from.
        If dict, {key:url} format and keys will be used in output dict.
        Of list, numerical keys will be used in output.

    dataName : str, optional, default = None
        Subdir to use for data.
        If None will default to 'downloads'

    dataPath : str or Path object, optional, default = None
        Full dir path to use.
        If None will default to current working dir.
        Note that 'dataName' is added to this path.

    ghRaw : bool, optional, default = True
        If true add '?raw=true' if 'github.com' is in the URL.

    Returns
    -------
    list : list of local files

    dict : dictionary or local files, including source URL and if freshly downloaded


    Notes
    -----

    01/05/23    v1  Fairly basic. Should add some optional flags for paths etc.

    """

    if dataName is None:
        dataName = 'downloads'

    # Set data dir
    if dataPath is None:
        dataPath = Path(Path.cwd(), dataName)

    if not isinstance(dataPath,Path):
        dataPath = Path(dataPath)

    # Create and pull files if dir not present (NOTE doesn't check per file here)
    if not dataPath.is_dir():
        dataPath.mkdir()
        print(f'*** Created output dir {dataPath}')

    # Set dict if list is passed
    if isinstance(urls,list):
        urls = {n:item for n,item in enumerate(urls)}

    # Check files and pull if missing

    fList = []
    fDict = {}

    for k,v in urls.items():

        # Check for local file
        localFile = dataPath/Path(v).parts[-1]

        if localFile.exists():
            print(f'Local file {localFile} already exists')

            fList.append(localFile)
            fDict[k] = {'url':v, 'localFile':localFile, 'downloaded':False}

        else:
            print(f'Downloading from {v} to {localFile}.')

            # Pull files with wget
            if ghRaw and (Path(v).parts[1] == 'github.com'):
                v = v+'?raw=true'    # For Github add '?raw=true' to URL

            fout = wget.download(v,out=dataPath.as_posix())

            # fout = wget.download(item['links']['self'], out=localFile.as_posix())

            print(f"Pulled to file: {fout}")

            fList.append(Path(fout))  # Log local file list

            fDict[k] = {'url':v, 'localFile':Path(fout), 'downloaded':True}

        # Add file info
        fDict[k]['filestat'] = fDict[k]['localFile'].lstat()
        fDict[k]['size'] = convert_bytes(fDict[k]['filestat'].st_size)
        fDict[k]['timestamp'] = dt.fromtimestamp(fDict[k]['filestat'].st_ctime).strftime('%H:%M:%S %d/%m/%Y')


    return fList, fDict
