# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:light
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.4'
#       jupytext_version: 1.2.4
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# # Set plot options - interactive to json file
# 07/08/20
#
# For wavefunction plotting, now changing format to use dict for plot options, and json file to store/modify. For ease of use & reproducibility (and to save writing the raw file), the original settings are set here.
#
# DEPRECATED - see end of notebook for current packaged function version.

# +
import json
import pprint
pp = pprint.PrettyPrinter(indent=4)

import sys
from pathlib import Path

# ePSproc test codebase (local)
# For package version this shouldn't be necessary
if sys.platform == "win32":
    modPath = r'D:\code\github\ePSproc'  # Win test machine
else:
    modPath = r'/home/femtolab/github/ePSproc/'  # Linux test machine

optionsFile = Path(modPath,'epsproc','vol','plotOptions.json')
# -


# ## Read existing options (if set)

optionsFile.is_file()

# +
if optionsFile.is_file():
    with open(optionsFile) as json_file:
        optionsFileJSON = json.load(json_file)
    
    print(f"*** Read existing file {optionsFile} OK, contents:")
    print(json.dumps(optionsFileJSON, sort_keys=False, indent=4))
    
else:
    print(f"*** File {optionsFile} not found.")
    

# -

# ## Set new options

optionsLocal = {}

# ### Define global settings

# +
globalSettings = {"note":"Global plot settings, used as defaults. To change for session, overwrite in local dict. To change permanently, overwrite in file plotOptions.json. To reset, use `epsproc/vol/set_plot_options_json.ipynb` or .py.",
                  "pType":"Abs", "interactive":True, "inline":True, "animate":False,
                  "isoLevels":6, "isoValsAbs":None, "isoValsPC":None,
                  "opacity":0.5,
#                   "plotter":""  # Set plotter dynamically based on options above...?
                  }

optionsLocal["global"] = globalSettings
# -

# ## Define plotter specific settings

# +
BGplotterSettings = {"addAxis" : True,
                     "kwargs" : {}    # Set empty kwargs dict for passing any other params at run time.
                    }

optionsLocal["BGplotter"] = BGplotterSettings
# -

# ## Write settings to file

with open(optionsFile, 'w') as json_file:
    
    if optionsFile.is_file():
        owFlag = input(f"File {optionsFile} exists, overwrite (y/n)? ")
    else:
        owFlag = 'y'
        
    if owFlag == 'y':
        json.dump(optionsLocal, json_file, indent=4, sort_keys=False)  # Set indent + sort keys for nicer (HF) file output.

stop here if running full notebook!!!

# ## Functional version testing...

# +
sys.path.append(modPath)
# import epsproc as ep

from epsproc.vol.setOptions import setLocalOptions, readOptionsFile, writeOptionsFile  # Plot options IO
# -

output = setLocalOptions()

output

writeOptionsFile(Path("test.out"), output)


