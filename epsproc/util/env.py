"""
ePSproc util functions for environment checks and related.

28/09/20    v1  Notebook checking fn. added, use for some types of display settings.

"""

# To control display of warnings...
import warnings


# Check if iPython/Notebook
# Code from: https://exceptionshub.com/how-can-i-check-if-code-is-executed-in-the-ipython-notebook.html
# FAILS - item doesn't exist any more it seems!
# def in_ipynb():
#     try:
#         cfg = get_ipython().config
#         if cfg['IPKernelApp']['parent_appname'] == 'ipython-notebook':
#             return True
#         else:
#             return False
#     except NameError:
#         return False

# THIS works, IPython v7.1.3, Jupyter notebook 6.0.3
# COULD also use Scooby here, see https://pypi.org/project/scooby/
def isnotebook(warn = 'once'):
    """
    Check if code is running in Jupyter Notebook.

    Taken verbatim from https://exceptionshub.com/how-can-i-check-if-code-is-executed-in-the-ipython-notebook.html
    Might be a better/more robust way to do this?

    Also added optional "warn = 'once'", set 'once' or 'ignore' to control warnings in notebook env., see https://stackoverflow.com/a/9031848

    """

    try:
        shell = get_ipython().__class__.__name__

        if shell == 'ZMQInteractiveShell':
            warnings.filterwarnings(warn)
            return True   # Jupyter notebook or qtconsole

        elif shell == 'TerminalInteractiveShell':
            return False  # Terminal running IPython

        else:
            return False  # Other type (?)

    except NameError:
        return False      # Probably standard Python interpreter
