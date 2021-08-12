###########################
Extended installation notes
###########################
05/11/20

ePSproc is available on PyPi, which is the simplest means of installing the latest stable release. For testing and/or development updates, use the Github repo source.

It has been tested with Python >= 3.7, and various package versions. Currently (v1.3.0, 12th August 2021), up-to-date packages are working apart from Xarray (>=0.12, <=0.15.0) and Seaborn (==0.9.0) - this will be fixed in future. (Note that packages in `requirements.txt` are currently unversioned apart from these, see notes below for more details.)

In all cases, you might want to use a virtual environment, e.g.:

With Anaconda

.. code-block::

  conda create --name epsProc python=3.7
  conda activate epsProc


With venv (see https://docs.python.org/3/library/venv.html)

.. code-block::

  python -m venv /path/to/new/virtual/environment


From PyPi
=========

Pip only
--------

The simplest method, this will pull the latest packaged release, and dependencies. (There is currently no Conda packaged version.)

``pip install ePSproc``

Note that some versions may not pull dependencies (not set in ``setup.py`` in some cases, should be OK for versions `ePSproc>=1.3.0-dev`), these can be additionally installed from the `requirements.txt file <https://github.com/phockett/ePSproc/blob/master/requirements.txt>`_:

``pip install -r requirements.txt``


Or manually... adding a couple of core packages should cover most dependencies:

.. code-block::

  pip install xarray
  pip install --user git+git://github.com/moble/spherical_functions


Conda
-----

Either as pip case above, plus Conda for the dependencies; or from source - see below.




From source + local pip install
===============================

Using pip
---------

.. code-block::

  git clone https://github.com/phockett/ePSproc.git
  pip install -e epsproc

This should install all dependencies (for `ePSproc>=1.3.0-dev`), although ``spherical_functions`` are currently giving issues with pip+PyPi (tested 12th August 2021), so may need to be installed separately from source with ``pip install git+git://github.com/moble/spherical_functions`` (this is `Moble's spherical functions library <https://github.com/moble/spherical_functions>`_, but will be updated to the newer `spherical <https://github.com/phockett/ePSproc/issues/35>`_ package).

To install with specific dependencies, just run ``pip install -e epsproc -r ePSproc/requirements.txt``


From `setup.py`

.. code-block::

  git clone https://github.com/phockett/ePSproc.git
  cd ePSproc
  python setup.py install


Note this currently installs without dependencies  for `ePSproc<1.3.0-dev`.



Notes

* For a single branch use ``git clone --single-branch --branch <branchname> https://github.com/phockett/ePSproc.git``
* The repo can be passed directly to pip, e.g. ``pip install git+https://github.com/phockett/ePSproc.git``, see `notes in the pip docs <https://pip.pypa.io/en/stable/reference/pip_install/#git>`_.
* Note that ``pip -e`` is for 'editable', and requires the source dir to remain, but the installation is also editable, `see notes here <https://stackoverflow.com/questions/41535915/python-pip-install-from-local-dir>`_. Drop the ``-e`` for a normal installation.
* In v1.3.0 (12/08/21), setup.py does not contain a pkg requirements list, but this will change in future (see `discussion here <https://stackoverflow.com/a/33685899>`_.)


Assuming a fresh environment, you might also need to install some requirements manually:

With pip
  pip install git+git://github.com/moble/quaternion
  pip install git+git://github.com/moble/spherical_functions


Using Conda
-----------

.. code-block::

  git clone https://github.com/phockett/ePSproc.git
  conda create --name ePSproc --file ePSproc/requirements.txt --channel default --channel conda-forge
  pip install -e epsproc


However... this may fail if any of the packages are missing or give issues. A quick fix is to `iterate over lines <https://stackoverflow.com/questions/35802939/install-only-available-packages-using-conda-install-yes-file-requirements-t>`_

.. code-block::

  conda create --name epsProc python=3.7
  conda activate epsProc
  while read requirement; do conda install --yes $requirement --channel default --channel conda-forge; done < ePSproc/requirements.txt
  pip install -e epsproc

Note that the python version is optional here, and the latest version will be pulled on install if not specified.


With specific (working) package versions: select from the various `*.yml` files under `/notes/envs`. These correspond to tested working environments.




Other options
===================

A few other options...

- Some development envs are available as Conda .yml files in ``/notes/envs``, these can be used to clone a known-working env.
  - E.g. ``conda env create -f environment_epsdev_v1.3.0_040821_no-builds.yml`` for the current v1.3.0 environment.
  - Note that these envs may included extra packages and/or platform specific packages.
  - Note that the envs are currently a bit of a mess, but will be cleared up soon.
  - See the `Conda docs for more details on .yml env files <https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#sharing-an-environment>`_

- Similarly, there are some alternative Pip requirement.txt files in ``/notes/envs``, corresponding to some specific development envs.
  - E.g. ``pip install -r requirements_epsdev_v1.3.0_040821.txt`` for the current v1.3.0 environment.
  - Note that these envs may included extra packages and/or platform specific packages.
  - Note that the envs are currently a bit of a mess, but will be cleared up soon.

- If using nb_conda_kernels:
  conda install ipykernel




--------------

For Conda:

  conda install -c conda-forge spherical_functions






RuntimeError: The current Numpy installation ('C:\\Users\\femtolab\\AppData\\Local\\Temp\\pip-build-env-gqcuairg\\overlay\\Lib\\site-packages\\numpy\\__init__.py') fails to pass a sanity check due to a bug in the windows runtime. See this issue for more information: https://tinyurl.com/y3dm3h86

May need first:
  conda install numpy numba

CURRENTLY issues with numpy v1.19.2 (Win), see https://developercommunity.visualstudio.com/content/problem/1207405/fmod-after-an-update-to-windows-2004-is-causing-a.html

* conda install numpy==1.19.1 didn't fix.
* conda install numpy==1.17.4 (random old version!) didn't fix.
* pip install numpy==1.19.3 didn't fix

Q: is pip pulling latest version at pip install?
