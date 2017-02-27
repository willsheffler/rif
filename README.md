tell willsheffler@gmail.com to put something in this readme!

python requirements
numpy, jinja2, pytest, pytest-cpp*, hypothesis, pandas
numpy-quaternion

*use the pytest\_cpp in external: "cd \<rifdi\r>/external/hacked_packages/pytest-cpp-0.4 && pip install ."

for graphics requirements
pymol + pip packages: PyOpenGL OpenGLContext PyVRML97 pydispatcher

SUBLIME SETUP
I use sublime plugins EasyClangComplete, Git, C++11, GitGutter, SublimeLinter

NOTES:
iff you have pyrosetta4 on your path, functions using pyrosetta will be used/tested
iff you have pymol on your path, functions using pymol will be used/tested

-------------- ubuntu 16.04 --------------
\# this ppa is to install gcc-6
sudo add-apt-repository ppa:ubuntu-toolchain-r/test # optional, for gcc6
sudo apt update
sudo sudo apt install git
git clone git@github.com:willsheffler/rif.git
sudo apt-get install python-pip python3-pip cmake ninja-build clang g++-6 libboost-system-dev libboost-iostreams-dev
sudo apt-get install clang g++-6 # optional
cd external/hacked_packages/pytest-cpp-0.4 && sudo -H pip2 install . && cd -
sudo -H pip2 install future pytest-xdist hypothesis jinja2 tox numpy pandas
sudo -H pip2 install numpy-quaternion
cd external/hacked_packages/pytest-cpp-0.4 && sudo -H pip3 install . && cd -
sudo -H pip3 install future pytest-xdist hypothesis jinja2 tox numpy pandas
sudo -H pip3 install numpy-quaternion
python tools/cmake_build_and_run_pytest.py # convenient script to build and test without installing
tox

----------------- mac 10.10 -------------------
brew install python cmake ninja boost
/usr/local/bin/pip2 install future pytest-xdist hypothesis jinja2 tox numpy pandas
/usr/local/bin/pip2 install numpy-quaternion
cd external/hacked_packages/pytest-cpp-0.4 && /usr/local/bin/pip2 install . && cd -
/usr/local/bin/python tools/cmake_build_and_run_pytest.py