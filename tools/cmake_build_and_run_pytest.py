#!/usr/bin/env python

from __future__ import print_function

import sys
from build_utils import build_and_run_pytest

if __name__ == '__main__':
    try:
        sys.exit(build_and_run_pytest(redo_cmake=True))
    except Exception as e:
        print('!!!!! ' * 20)
        print(e)
        print('build_and_run_pytest exception! (above)')
        print('!!!!! ' * 20)
        sys.exit(-1)
