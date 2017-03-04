#!/usr/bin/env python

import sys
from build_utils import build_and_run_pytest

if __name__ == '__main__':
    sys.exit(build_and_run_pytest(redo_cmake=True))
