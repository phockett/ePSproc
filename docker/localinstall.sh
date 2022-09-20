#!/bin/sh

# Quick local install for ePSproc and PEMtk
# NOTE: this is working for root user, but pemtk has permissions issue for jovyan --user case, not sure why.

BASEPATH="${1:-/home/jovyan/github}"
# USER="${2:}"   # TODO: set user flag as option

cd $BASEPATH
pip install -e epsproc #--user
pip install -e pemtk #--user
