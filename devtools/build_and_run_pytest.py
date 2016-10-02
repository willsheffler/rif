from __future__ import absolute_import, division, print_function
from builtins import *

import os
import sys

os.sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from build_utils import get_proj_root, get_build_dir, rebuild_fast, add_to_pypath
from build_utils import remove_installed_riflib

if __name__ == '__main__':
    # remove_installed_rif()
    proj_root = get_proj_root()
    print('calling rebuild_fast')
    if rebuild_fast(target='riflib gtest_all', cfg='Release'):
        sys.exit(-1)
    pypath = [
        os.path.abspath(get_build_dir('lib')),
        proj_root + '/external',
    ]
    # need to use sys.path for this process
    sys.path = pypath + sys.path
    # need to use PYTHONPATH env for xdist subprocesses
    add_to_pypath(pypath)
    # TODO both here and in docs, this gets messed up when riflib is actually installed
    import pytest
    if sys.version_info.major is 2:
        proj_root = bytes(proj_root, 'ascii')
    pytest.main(proj_root)

