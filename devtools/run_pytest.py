import os
import sys
import pytest

os.sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from build_utils import get_proj_root, get_build_dir, add_to_pypath


if __name__ == '__main__':
    proj_root = get_proj_root()
    pypath = [
        os.path.abspath(get_build_dir('lib')),
        proj_root + '/external',
        ]
    # need to use sys.path for this process
    sys.path = pypath + sys.path
    # need to use PYTHONPATH env for xdist subprocesses
    add_to_pypath(pypath)

    pytest.main(proj_root)
