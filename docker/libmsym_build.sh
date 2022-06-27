#!/bin/sh

# Basic build script for libmsym, adapted from https://github.com/mcodev31/libmsym#installing
# Requires Cmake

git clone https://github.com/mcodev31/libmsym.git
cd libmsym
mkdir build
cd build
# build as shared library; build examples (built in ./examples,  not installed)
cmake -DBUILD_SHARED_LIBS:BOOL=ON -DMSYM_BUILD_EXAMPLES:BOOL=ON ../.
make
# sudo only required if installing in directory not owned by user
# use -DCMAKE_INSTALL_PREFIX:PATH=<libmsym installation path> to change

# May need sudo
make install

# Python install
cd ../bindings/python
# install libmsym module in user site
python setup.py install


# Set lib path
ldconfig /usr/local/lib/
