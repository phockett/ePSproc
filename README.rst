.. Readme originally converted from readme.md, via Pandoc
   pandoc -s -o readme.rst README.md

ePSproc Readme
==============

|docs|

A post-processing suite for ePolyScat calculations, and general tool for photoionzation calculations.

ePSproc scripts are designed for photoionization studies. The scripts were originally written for Matlab (2009 - 2016); a Python version is currently under (heavy) development (2019 to date), and is now the main/preferred version.

Source code is `available on Github <https://github.com/phockett/ePSproc>`_.

Ongoing documentation is on `Read the Docs <https://epsproc.readthedocs.io>`_.

For more background, and details on the Matlab version, see the software metapaper for ePSproc (Aug. 2016), *ePSproc: Post-processing suite for ePolyScat electron-molecule scattering calculations*, on `Authorea <https://www.authorea.com/users/71114/articles/122402/_show_article>`_ or `arXiv 1611.04043 <https://arxiv.org/abs/1611.04043>`_. For data platform and related experimental tools, see `the Photoelectron Metrolgy toolkit <https://pemtk.readthedocs.io/en/latest/about.html>`__.

.. image:: https://epsproc.readthedocs.io/en/latest/_images/output_12_2.png


Installation (Python)
---------------------

Basics:

* From source: simply download `from Github <https://github.com/phockett/ePSproc>`_ and use locally, or run ``pip install git+https://github.com/phockett/ePSproc.git``.

* From Pypi ``pip install ePSproc``

See the `extended installation notes for more details & options <etc/installation_notes_051120.html>`_.


Python
------

Functionality:

* Plot basic properties (cross-secions, anisotropy parameters) from ePS ``GetCro`` outputs.
* Read raw photoionization matrix elements from ePS output files with ``dumpIdy`` segments.
* Process matrix elements:

  * Plot matrix elements.
  * Calculate MF-PADs from the matrix elements.
  * Calculate MF, LF and AF $\beta_{LM}$ from the matrix elements.
  * Plot $\beta_{LM}$ & PADs.

* `Distirbution via PyPi (latest stable version) <https://pypi.org/project/ePSproc/>`__ .`
* `Docs on Read the Docs <https://epsproc.readthedocs.io>`__ .

.. This doesn't work for PyPi: See the demo :doc:`Jupyter notebook <ePSproc_demo_Aug2019/ePSproc_demo_Aug2019>` for example usage.

To get started, see the demo Jupyter notebooks for example usage:

* `Base class usage demo <https://epsproc.readthedocs.io/en/latest/demos/ePSproc_class_demo_161020.html>`__ (note this wraps most functionality as of Oct. 2020).
* `Basic functional demo <https://epsproc.readthedocs.io/en/latest/ePSproc_demo_Aug2019/ePSproc_demo_Aug2019.html>`__ .
* `MF $\beta_{LM}$ parameters demo <https://epsproc.readthedocs.io/en/latest/ePSproc_BLM_calc_demo_Sept2019_rst/ePSproc_BLM_calc_demo_Sept2019.html>`__ .



Source:

* ./epsproc: python version, code still under development.

* ./docs: documentation tree, `HTML version on Read the Docs <https://epsproc.readthedocs.io>`__.

.. This doesn't work for PyPi :doc:`Full function documentation <modules/epsproc>`.

`Full function documentation <https://epsproc.readthedocs.io/en/latest/modules/epsproc.html>`_.


Matlab
------

Functionality:

* Read raw photoionization matrix elements from ePS output files with ``dumpIdy`` segments
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

For more details see

* `ePSproc Matlab demo notebook <https://epsproc.readthedocs.io/en/latest/demos/ePSproc_Matlab_demo_notebook_090821.html>`_ .
* `ePSproc: Post-processing suite for ePolyScat electron-molecule scattering calculations <https://www.authorea.com/users/71114/articles/122402/_show_article>`_ .


Resources
---------

An ongoing repository of `ePS results can be found at ePSdata <https://phockett.github.io/ePSdata/index.html>`_ (as of Jan 2020, this replaces `the previous repository on OSF <https://osf.io/psjxt/>`_).



ePolyScat
---------

For details about ePolyScat (ePS), a tool for computation of electron-molecule scattering, see:

* `ePS website & manual <https://epolyscat.droppages.com>`_, maintained by R.R. Lucchese.

* Calculation of low-energy elastic cross sections for electron-CF4 scattering, F. A. Gianturco, R. R. Lucchese, and N. Sanna, J. Chem. Phys. 100, 6464 (1994), http://dx.doi.org/10.1063/1.467237

* Cross section and asymmetry parameter calculation for sulfur 1s photoionization of SF6, A. P. P. Natalense and R. R. Lucchese, J. Chem. Phys. 111, 5344 (1999), http://dx.doi.org/10.1063/1.479794


Future aims
-----------

-  Add capabilities, including more general processing, and for other phenomena (e.g. recombination matrix elements for high-harmonic generation, aligned-frame calculations)
-  Tidy and streamline code (yep)
-  Extend & update notes and benchmark calculations
-  Port to non-commercial run-time engines, e.g. python

Citation
--------

If you make use of ePSproc in your research, please cite it.

Cite the software directly via either Github or Figshare repositories for the software (note same DOI for both)::

  @misc{ePSprocGithub,
    title={ePSproc: Post-processing for ePolyScat},
    url={https://github.com/phockett/ePSproc},
    DOI={10.6084/m9.figshare.3545639},
    publisher={Github},
    howpublished = {\url{https://github.com/phockett/ePSproc}},
    author={Hockett, Paul},
    year={2016},
    commit = {30158eb3fbba41d0a4c3a973744f28b7187e6ee2}
  }

  @misc{ePSprocFigshare,
    title={ePSproc: Post-processing for ePolyScat},
    url={https://figshare.com/articles/ePSproc_Post-processing_for_ePolyScat_v1_0_0_/3545639/4},
    DOI={10.6084/m9.figshare.3545639},
    publisher={Figshare},
    author={Hockett, Paul},
    year={2016}
  }

... or the software paper (Authorea/arXiv)::

  @article{ePSprocPaper,
    title={ePSproc: Post-processing for ePolyScat electron-molecule scattering calculations},
    url={https://www.authorea.com/users/71114/articles/122402-epsproc-post-processing-suite-for-epolyscat-electron-molecule-scattering-calculations},
    DOI={10.22541/au.156754490.06103020},
    journal = {Authorea/arXiv e-prints},
    publisher={Authorea/arXiv},
    author={Hockett, Paul},
    year={2016},
    archivePrefix = {arXiv},
    eprint = {1611.04043},
    primaryClass = {physics.comp-ph},
    eid = {arXiv:1611.04043},
    pages = {arXiv:1611.04043}
  }

(Citation styles for software `from StackExchange <https://academia.stackexchange.com/questions/14010/how-do-you-cite-a-github-repository>`_.)

.. .. include:: citation.txt (keep duplicate details here, since this doesn't work for basic Github readme!)

Acknowledgements
----------------

Special thanks to R.R. Lucchese and coworkers for `ePolyScat <https://epolyscat.droppages.com>`_.

Thanks to the multiple collaborators and co-authors who encouraged and suggested the cavilier use of ePS "out of the box", for many different problems incorporating electron scattering and photoionization. This spirit of "shoot first, ask questions later" indeed raised many questions which, eventually, necessitated proper use of ePS and careful post-processing of the results, and sharpened related foundational expertise - efforts well worth making.

Thanks, finally, and of course, to those supporting scientific software development and infrastructure (and making it easy!), including Github, Read the Docs, Pypi, SciPy etc. etc. In particular the python version of this project makes use of `Xarray <http://xarray.pydata.org/en/stable/index.html>`_, and `Moble's spherical functions (& quaternion) <https://github.com/moble/spherical_functions>`_.

.. |docs| image:: https://readthedocs.org/projects/epsproc/badge/?version=latest
    :alt: Documentation Status
    :scale: 100%
    :target: https://epsproc.readthedocs.io/en/latest/?badge=latest
