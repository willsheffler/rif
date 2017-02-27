tell willsheffler@gmail.com to put something in this readme!

python requirements
numpy, jinja2, pytest, pytest-cpp*, hypothesis, pandas
numpy-quaternion

*use the pytest\_cpp in external: "cd \<rifdi\r>/external/hacked_packages/pytest-cpp-0.4 && pip install ."

for graphics requirements
pymol + pip packages: PyOpenGL OpenGLContext PyVRML97 pydispatcher

SUBLIME SETUP
I use sublime plugins EasyClangComplete, Git, C++11, GitGutter, SublimeLinter

-------------- ubuntu 16.04 --------------
sudo add-apt-repository ppa:ubuntu-toolchain-r/test
sudo apt update
sudo sudo apt install git
git clone git@github.com:willsheffler/rif.git
sudo apt-get install python-pip cmake ninja-build clang gcc-6 libboost-system-dev libboost-iostreams
sudo -H pip2 install future pytest pytest-xdist pytest-sugar hypothesis
