from __future__ import absolute_import, division, print_function, unicode_literals
from builtins import *

import os
import sys

os.sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from build_utils import get_proj_root, get_build_dir, rebuild_fast, add_to_pypath
from build_utils import remove_installed_riflib


if __name__ == '__main__':
    # remove_installed_rif()
    proj_root = get_proj_root()
    if rebuild_fast(target='gtest_wip'):
        sys.exit(-1)
    os.chdir(os.path.abspath(get_build_dir('temp')))
    os.system('./gtest_wip')


