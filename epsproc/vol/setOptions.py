
"""
ePSproc vol module setOptions

Functions to read & write default set of plotting options to file.

If run as main:

- Check existing file from passed arg, or in default location (epsproc/vol/plotOptions.json)
- Read file if exists.
- If file is missing, prompt to write defaults to file.

08/08/20    v1, dev.  See also set_plot_options_json.ipynb

"""

import json
import pprint
pp = pprint.PrettyPrinter(indent=4)

import sys
import os
import inspect
from pathlib import Path



def setLocalOptions():
    optionsLocal = {}

    globalSettings = {"note":"Global plot settings, used as defaults. To change for session, overwrite in local dict. To change permanently, overwrite in file plotOptions.json. To reset, use `epsproc/vol/set_plot_options_json.ipynb` or .py.",
                      "pType":"Abs", "interactive":True, "inline":True, "animate":False,
                      "isoLevels":6, "isoValsAbs":None, "isoValsPC":None, "isoValsGlobal":True,
                      "opacity":0.5,
    #                   "plotter":""  # Set plotter dynamically based on options above...?
                      }

    optionsLocal["global"] = globalSettings


    BGplotterSettings = {"addAxis" : True,
                         "kwargs" : {}    # Set empty kwargs dict for passing any other params at run time.
                        }

    optionsLocal["BGplotter"] = BGplotterSettings

    return optionsLocal


# def setOptionsFile(optionsFile = None):

def readOptionsFile(optionsFile, verbose = False):

    # Set path wrapper in case str was passed.
    optionsFile = Path(optionsFile)

    if optionsFile.is_file():
        # try:
        with open(optionsFile) as json_file:
            optionsFileJSON = json.load(json_file)

        print(f"\n*** Read existing plot options from file {optionsFile} OK.")
        if verbose:
            print(json.dumps(optionsFileJSON, sort_keys=False, indent=4))

        return optionsFileJSON

    else:
        print(f"\n*** Plot options file {optionsFile} not found, using defaults.")

        return setLocalOptions()


def writeOptionsFile(optionsFile, optionsLocal, owFlag = False):
    # Set path wrapper in case str was passed.
    optionsFile = Path(optionsFile)

    print(f"*** Writing plot options to file {optionsFile}")

    if optionsFile.is_file():
        owFlag = input(f"File {optionsFile} exists, overwrite (y/n)? ")
    else:
        owFlag = 'y'

    with open(optionsFile, 'w') as json_file:
        if owFlag == 'y':
            json.dump(optionsLocal, json_file, indent=4, sort_keys=False)  # Set indent + sort keys for nicer (HF) file output.


if __name__ == "__main__":

    # Check passed args
    if len(sys.argv > 1):
        optionsFile = Path(sys.argv[1])
    else:
        optionsFile = None

    # Set default path based on file location - may not be robust?
    # From https://stackoverflow.com/a/12154601
    if optionsFile is None:
        optionsFile = Path((os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))),'plotOptions.json')
        writeFlag = True

    # Read file
    optionsLocal = setLocalOptions()

    # if optionsFile.is_file():
    #     with open(optionsFile) as json_file:
    #         optionsFileJSON = json.load(json_file)
    #
    #     print(f"*** Read existing file {optionsFile} OK, contents:")
    #     print(json.dumps(optionsFileJSON, sort_keys=False, indent=4))
    #
    # else:
    #     print(f"*** File {optionsFile} not found.")

    if writeFlag:
        ow = input("Write defaults to {optionsFile} (y/n)? ")

        if ow == 'y':
            writeOptionsFile(optionsFile, setLocalOptions())
