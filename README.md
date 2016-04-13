# ePSproc
Post-processing suite for ePolyScat calculations.

ePSproc scripts are designed for photoionization studies, and specifically:
- Read raw photoionization matrix elements from ePS output files with "dumpIdy" segments
- Calculate MF-PADs from the matrix elements
- Plot MF-PADs
- Plot X-sects

The scripts are currently written for Matlab. The distribution currently contains a set of functions for processing, a script showing demo calculations, and the benchmark results from these calculations.

# ePolyScat
For details about ePolyScat (ePS), a tool for computation of electron-molecule scattering, see:
- ePS website & manual, maintained by R.R. Lucchese, http://www.chem.tamu.edu/rgroup/lucchese/ePolyScat.E3.manual/manual.html
- F. A. Gianturco, R. R. Lucchese, and N. Sanna, J. Chem. Phys. 100, 6464 (1994), http://dx.doi.org/10.1063/1.467237
- A. P. P. Natalense and R. R. Lucchese, J. Chem. Phys. 111, 5344 (1999), http://dx.doi.org/10.1063/1.479794

# Future aims
- Add capabilities, including more general processing, and for other phenomena (e.g. recombination matrix elements for high-harmonic generation, aligned-frame calculations)
- Tidy and streamline code (yep)
- Port to non-commercial run-time engines, e.g. python

