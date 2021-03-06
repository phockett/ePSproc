# ePSproc
Post-processing suite for ePolyScat calculations.

PKG branch - testing package structure (Hatch), docs (Sphinx) and modules.

ePSproc scripts are designed for photoionization studies, and specifically:
- Read raw photoionization matrix elements from ePS output files with "dumpIdy" segments
- Calculate MF-PADs from the matrix elements (ePSproc_MFPAD.m, see also ePSproc_NO2_MFPADs_demo.m)
- Plot MF-PADs
- Plot X-sects
- (Beta testing): Calculate MF-BLMs from matrix elements, see ePSproc_MFBLM.m
- (Under development): Calculate AF-BLMs from matrix elements.

The scripts are currently written for Matlab & Python. The distribution currently contains:
- ./matlab: stable matlab code (as per [release v1.0.1](https://github.com/phockett/ePSproc/releases)).
- ./python-dev: basic python version, code still under development.

The Matlab distro contains:
- a set of functions for processing (ePSproc*.m files)
- a script showing demo calculations, <a href="https://github.com/phockett/ePSproc/blob/master/ePSproc_NO2_MFPADs_demo.m">ePSproc_NO2_MFPADs_demo.m</a>, and the benchmark results from these calculations, <a href="https://github.com/phockett/ePSproc/blob/master/ePSproc_NO2_testing_summary_250915.pdf">ePSproc_NO2_testing_summary_250915.pdf</a>
- some additional notes on ePS photoionization matrix elements, <a href="https://github.com/phockett/ePSproc/blob/master/ePSproc_scattering_theory_ePS_notes_011015.pdf">ePSproc_scattering_theory_ePS_notes_011015.pdf</a>.

A software metapaper for ePSproc can be found on Authorea: https://www.authorea.com/users/71114/articles/122402/_show_article

An ongoing repository of ePS results can be found on OSF: https://osf.io/psjxt/

# ePolyScat
For details about ePolyScat (ePS), a tool for computation of electron-molecule scattering, see:
- ePS website & manual, maintained by R.R. Lucchese, https://epolyscat.droppages.com
- F. A. Gianturco, R. R. Lucchese, and N. Sanna, J. Chem. Phys. 100, 6464 (1994), http://dx.doi.org/10.1063/1.467237
- A. P. P. Natalense and R. R. Lucchese, J. Chem. Phys. 111, 5344 (1999), http://dx.doi.org/10.1063/1.479794

# Future aims
- Add capabilities, including more general processing, and for other phenomena (e.g. recombination matrix elements for high-harmonic generation, aligned-frame calculations)
- Tidy and streamline code (yep)
- Extend & update notes and benchmark calculations
- Port to non-commercial run-time engines, e.g. python
