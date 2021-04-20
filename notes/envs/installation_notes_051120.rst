
Extended installation notes
05/11/20

ePSproc is available on PyPi, which is the simplest means of installing the latest stable release. For testing and/or development updates, use the Github repo source.

In all cases, you might want to use a virtual environment, e.g.:

With Anaconda

.. code-block::

  conda create --name epsProc --python=3.7
  conda activate epsProc


With venv (see https://docs.python.org/3/library/venv.html)

.. code-block::

  python -m venv /path/to/new/virtual/environment


From source + local pip install

Using pip

.. code-block::

  git clone https://github.com/phockett/ePSproc.git
  pip install -e epsproc

[No dependencies installed in this case]

From setup.py

``python setup.py install``



Notes

* For a single branch use `git clone --single-branch --branch <branchname> https://github.com/phockett/ePSproc.git`
* The repo can be passed directly to pip, e.g. `pip install git+https://github.com/phockett/ePSproc.git`, see `notes in the pip docs <https://pip.pypa.io/en/stable/reference/pip_install/#git>`_.
* Note that `pip -e` is for 'editable', and requires the source dir to remain, but the installation is also editable, `see notes here <https://stackoverflow.com/questions/41535915/python-pip-install-from-local-dir>`_. Alternatively, use `pip install <path_to_local_pkg>`.


Assuming a fresh environment, you might also need to install some requirements manually:

For Conda:

conda install -c moble spherical_functions

https://github.com/moble/spherical_functions
CURRENTLY missing channel??? Conda won't find this or numpy_quaternion

With pip
  pip install git+git://github.com/moble/quaternion
  pip install git+git://github.com/moble/spherical_functions

RuntimeError: The current Numpy installation ('C:\\Users\\femtolab\\AppData\\Local\\Temp\\pip-build-env-gqcuairg\\overlay\\Lib\\site-packages\\numpy\\__init__.py') fails to pass a sanity check due to a bug in the windows runtime. See this issue for more information: https://tinyurl.com/y3dm3h86

May need first:
  conda install numpy numba

CURRENTLY issues with numpy v1.19.2 (Win), see https://developercommunity.visualstudio.com/content/problem/1207405/fmod-after-an-update-to-windows-2004-is-causing-a.html

* conda install numpy==1.19.1 didn't fix.
* conda install numpy==1.17.4 (random old version!) didn't fix.
* pip install numpy==1.19.3 didn't fix

Q: is pip pulling latest version at pip install?
