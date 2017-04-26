#!/usr/bin/env python

from __future__ import print_function

import sys
import os
import traceback
sys.path.append(os.path.dirname(__file__))  # not sure why sometimes necessary
from build_utils import build_and_test, build_and_run_gtest_auto, _VERBOSE


if __name__ == '__main__':
    if _VERBOSE:
        print('== build_and_test.py start ==')
        print('== in ci ==')
        try:
            sys.exit(build_and_test())
        except Exception as e:
            t, v, tb = sys.exc_info()
            with open('.ERROR', 'w') as out:
                out.write(str(e))
            print('==========================================================')
            print("error running build_and_test, traceback:")
            print('==========================================================')
            print(e)
            traceback.print_tb(tb)
    else:
        try:
            build_and_run_gtest_auto()
            gtest_ran = True
        except NotImplementedError:
            gtest_ran = False
        # print("gtest_ran:", gtest_ran)
        if not gtest_ran:
            build_and_test()
