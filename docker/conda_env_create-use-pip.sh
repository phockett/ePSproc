#!/bin/bash

# Basic config for Conda env for ePSproc - VERSION WITH PIP INSTALLS
# 13/06/24
#
# Note this basically duplicates the steps in the current Docker builds, but conda/pip installs only.
# Up-to-date Docker build is QM3 version, https://github.com/phockett/Quantum-Metrology-with-Photoelectrons-Vol3/tree/main/docker
#
# NOTE: this does not currently build libmsym, see https://github.com/phockett/open-photoionization-docker-stacks/blob/main/epsproc-pemtk/libmsym_build.sh
#
# NOTE 2: this takes a LONG time, use pip install version for faster builds.
#

#*** Env name
ENVNAME="${1:-eps2024}"
BASEPATH="${2:-~/github}"   # Used only for local install option with existing local pkgs

#*** Set any versions here
PYTHON_VERSION=3.10.11    # In testing have some issues with 3.10.13, but this version is OK
CMAKE_VERSION=3.21.0
XARRAY_VERSION=2022.3.0   # Currently have IO issues for newer versions, see https://github.com/phockett/ePSproc/issues/64
SEABORN_VERSION=0.9       # Required for lmplot routine, see https://github.com/phockett/ePSproc/issues/23


#*** Create Conda env
echo Creating new conda env ${ENVNAME}

# From ePSproc requirements file
# FAILS... need to add channels? E.g. "conda install -f -y -q --name py33 -c conda-forge --file requirements.txt"
# Or can do line-by-line force, "while read requirement; do conda install --yes $requirement; done < requirements_minimal_unversioned.txt"
# conda env create -f {epsprocPath}/requirements.txt -n {ENVNAME}

# Bare env with installs below
conda create -y -n ${ENVNAME} python=${PYTHON_VERSION}

#*** Conda & pip installs
# Conda version
# conda install --quiet --yes -n ${ENVNAME} holoviews seaborn=${SEABORN_VERSION} selenium h5netcdf
# conda install --quiet --yes -n ${ENVNAME} -c conda-forge spherical_functions scooby jupytext pyshtools lmfit firefox geckodriver qutip xyzpy xarray=${XARRAY_VERSION}
# conda install --quiet --yes -n ${ENVNAME} -c pyviz hvplot

# Pip version - note assumes path here.
conda run -n ${ENVNAME} pip install -r ${BASEPATH}/ePSproc/requirements.txt

#*** Jupyterlab Plotly support
# Optional, run https://github.com/phockett/open-photoionization-docker-stacks/blob/main/epsproc-pemtk/plotlyinstall.sh
# plotlyinstall.sh
# conda install --quiet --yes -n ${ENVNAME} -c conda-forge -c plotly jupyter-dash


#*** Jupyter Book and build tools
# Optional
# pip install jupyter-book ghp-import jupyterlab-spellchecker
# pip install wget

# ****** ePSproc + PEMtk pkg installs

# *** Local PKG installs
# Optional, see https://github.com/phockett/open-photoionization-docker-stacks/blob/main/epsproc-pemtk/localinstall.sh

echo ***Installing local pkgs from ${BASEPATH}
conda run -n ${ENVNAME} sh ./localinstall.sh ${BASEPATH}

# *** PKG installs from GH
# Optional, see https://github.com/phockett/Quantum-Metrology-with-Photoelectrons-Vol3/blob/main/docker/scripts/ghinstall.sh
# ghinstall.sh


# *** Additional options and tools

# iPython Kernel for env
# conda install --quiet --yes -n ${ENVNAME} ipykernel
conda run -n ${ENVNAME} pip install ipykernel
