# ePSproc
Post-processing of ePolyScat calculations.

ePolyScat (ePS) is a tool for computation of electron-molecule scattering, distributed by R. R. Lucchese (http://www.chem.tamu.edu/rgroup/lucchese/ePolyScat.E3.manual/manual.html) .

ePSproc is a suite of scripts for post-processing of ePS results.  

ePSproc scripts are designed for photoionization studies, and specifically:
- Read raw photoionization matrix elements from ePS output files with "dumpIdy" segments
- Calculate MF-PADs from the matrix elements
- Plot MF-PADs
- Plot partial waves
- Plot X-sects

The scripts are currently written for Matlab.

# Future aims
- Add capabilities, including more general processing, and for other phenomena (e.g. high-harmonic generation recombination matrix elements)
- Tidy and streamline code (yep)
- Port to non-commercial run-time engines, e.g. python
