#!/usr/bin/env python

from __future__ import print_function

import sys
import os
from build_utils import build_and_test

if __name__ == '__main__':
    print('== build_and_test.py start ==')
    if 'CI' in os.environ:
        try:
            sys.exit(build_and_test())
        except Exception as e:
            print("error running build_and_test")
            with open('.ERROR', 'w') as out:
                out.write(str(e))
    else:
        build_and_test()
