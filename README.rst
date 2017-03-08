
Rotamer Interaction Field
=========================

.. inclusion-marker-do-not-remove


.. image:: https://img.shields.io/pypi/v/rif.svg
    :target: https://pypi.python.org/pypi/rif/

**master**

.. image:: https://img.shields.io/travis/willsheffler/rif.svg
    :target: http://travis-ci.org/willsheffler/rif
.. image:: https://img.shields.io/codecov/c/github/willsheffler/rif.svg
    :target: https://codecov.io/gh/willsheffler/rif


**devel**

.. image:: https://img.shields.io/travis/willsheffler/rif/devel.svg
    :target: http://travis-ci.org/willsheffler/rif
.. image:: https://img.shields.io/codecov/c/github/willsheffler/rif/devel.svg
    :target: https://codecov.io/gh/willsheffler/rif/devel


tell willsheffler@gmail.com to put something in this readme!
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


INSTALL
-------

python requirements
numpy, jinja2, pytest, pytest-cpp*, hypothesis, pandas

*use the pytest\_cpp in external: "cd \<rifdi\r>/external/hacked_packages/pytest-cpp-0.4 && pip install ."

for graphics requirements
pymol + pip packages: PyOpenGL OpenGLContext PyVRML97 pydispatcher

SUBLIME SETUP
~~~~~~~~~~~~~~
I use sublime plugins Anaconda, EasyClangComplete, Git, C++11, GitGutter, SublimeLinter

NOTES
~~~~~~~
iff you have pyrosetta4 on your path, functions using pyrosetta will be used/tested
iff you have pymol on your path, functions using pymol will be used/tested

docs
~~~~
install doxygen (--with-libclang probably best...)
pip3 install sphinx_rtd_theme breathe

ubuntu 16.04
~~~~~~~~~~~~
\# better to use a virtualenv, no sudo
\# this ppa is to install gcc-6
sudo add-apt-repository ppa:ubuntu-toolchain-r/test # optional, for gcc6
sudo apt update
sudo sudo apt install git
git clone git@github.com:willsheffler/rif.git
sudo apt-get install python-pip python3-pip cmake ninja-build clang g++-6 libboost-system-dev libboost-iostreams-dev
sudo apt-get install clang g++-6 # optional


mac (tested on 10.10 and 10.12)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
brew install python cmake ninja boost

universal
~~~~~~~~~~

pip install -rrequirements.txt
python tools/cmake_build_and_run_pytest.py

