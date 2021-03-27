.. Readme originally converted from readme.md, via Pandoc
   pandoc -s -o readme.rst README.md

ePSproc Readme
==============

Post-processing suite for ePolyScat calculations.

ePSproc scripts are designed for photoionization studies. The scripts were originally written for Matlab (2009 - 2016); a Python version is currently under (heavy) development (Aug. 2019), and will be the main version in future.

Source code is `available on Github <https://github.com/phockett/ePSproc>`_.

Ongoing documentation is on `Read the Docs <https://epsproc.readthedocs.io>`_.

For more background, and details on the Matlab version, see the software metapaper for ePSproc (Aug. 2016), *ePSproc: Post-processing suite for ePolyScat electron-molecule scattering calculations*, on `Authorea <https://www.authorea.com/users/71114/articles/122402/_show_article>`_ or `arXiv 1611.04043 <https://arxiv.org/abs/1611.04043>`_.

Installation
------------

From source: simply download `from Github <https://github.com/phockett/ePSproc>`_. See specific version notes below for more details on the source code.

Python:

.. code-block:: bash

    $ pip install ePSproc


Python
------

Functionality:

* Read raw photoionization matrix elements from ePS output files with “dumpIdy” segments
* Calculate MF-PADs from the matrix elements (ePSproc_MFPAD.m, see also ePSproc_NO2_MFPADs_demo.m)
* Plot MF-PADs
* Under development: additional functionality and distribution via PyPi.

See the demo :doc:`Jupyter notebook <ePSproc_demo_Aug2019/ePSproc_demo_Aug2019>` for example usage.

Source:

* ./epsproc: basic python version, code still under development.

* ./docs: documentation tree, `HTML version on Read the Docs <https://epsproc.readthedocs.io>`__.

:doc:`Full function documentation <modules/epsproc>`.


Matlab
------

Functionality:

* Read raw photoionization matrix elements from ePS output files with “dumpIdy” segments
* Calculate MF-PADs from the matrix elements (ePSproc_MFPAD.m, see also ePSproc_NO2_MFPADs_demo.m)
* Plot MF-PADs
* Plot X-sects
* (Beta testing): Calculate MF-BLMs from matrix elements, see ePSproc_MFBLM.m
* (Under development): Calculate AF-BLMs from matrix elements.


Source:

* /matlab: stable matlab code (as per `release v1.0.1 <https://github.com/phockett/ePSproc/releases>`__).

  * a set of functions for processing (ePSproc*.m files)
  * a script showing demo calculations, ePSproc_NO2_MFPADs_demo.m


* /docs/additional contains:

  * the benchmark results from these calculations, ePSproc_NO2_testing_summary_250915.pdf
  * additional notes on ePS photoionization matrix elements, ePSproc_scattering_theory_ePS_notes_011015.pdf.

See `ePSproc: Post-processing suite for ePolyScat electron-molecule scattering calculations <https://www.authorea.com/users/71114/articles/122402/_show_article>`_ for more details.


Resources
---------

An ongoing repository of `ePS results can be found on OSF <https://osf.io/psjxt/>`_.


ePolyScat
---------

For details about ePolyScat (ePS), a tool for computation of
electron-molecule scattering, see:

* ePS website & manual, maintained by R.R. Lucchese, https://epolyscat.droppages.com
* F. A. Gianturco, R. R. Lucchese, and N. Sanna, J. Chem. Phys. 100, 6464 (1994), http://dx.doi.org/10.1063/1.467237
* A. P. P. Natalense and R. R. Lucchese, J. Chem. Phys. 111, 5344 (1999), http://dx.doi.org/10.1063/1.479794


Future aims
-----------

-  Add capabilities, including more general processing, and for other phenomena (e.g. recombination matrix elements for high-harmonic generation, aligned-frame calculations)
-  Tidy and streamline code (yep)
-  Extend & update notes and benchmark calculations
-  Port to non-commercial run-time engines, e.g. python
