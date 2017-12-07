# gcc6, cmake
wget --no-check-certificate https://www.dropbox.com/s/bxc7ypidne8lw4i/centos5boost165.tar.gz?dl=0 -O - | tar -C / -xz
wget --no-check-certificate https://www.dropbox.com/s/15yqynnjj10k09c/cmake39centos5.tar.bz2?dl=0 -O - | tar -C / -xj
wget https://github.com/squeaky-pl/centos-devtools/releases/download/6.3/gcc-6.3.0-binutils-2.27-x86_64.tar.bz2 -O - | tar -C / -xj
export CC=/opt/devtools-6.3/bin/gcc
export CXX=/opt/devtools-6.3/bin/g++

git clone https://github.com/willsheffler/rif
cd rif
git checkout devel

export OLDPATH=$PATH

export PATH=/opt/python/cp35-cp35m/bin:/opt/devtools-6.3/bin:$OLDPATH
pip install -r requirements.txt
python setup.py bdist_wheel

export PATH=/opt/python/cp36-cp36m/bin:/opt/devtools-6.3/bin:$OLDPATH
pip install -r requirements.txt
python setup.py bdist_wheel


