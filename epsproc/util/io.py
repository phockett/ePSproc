"""
ePSproc IO util functions

Various tools for use in file IO.

27/06/22    Split out from core IO.py, to extend backend options and support.
            Now additionally wrapped therein for flexible handling of multiple backends.

"""

from datetime import datetime as dt # Import datetime.datetime for now() function

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
