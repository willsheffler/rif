#!/usr/bin/env python

"""build and run gtest_wip c++ exe"""

from __future__ import print_function
from builtins import *

import os
import sys

from build_utils import get_cmake_dir, rebuild_fast


def make_gtest_auto_cpp(files, cmake_dir):
    includes = "\n".join("#include <%s>" % f for f in files)
    code = """#include "test/gtest_util.hpp"
%s
int main(int argc, char **argv) {
  std::vector<std::string> args;
  for (int i = 0; i < argc; ++i) args.push_back(std::string(argv[i]));
  std::cout << "int main(int argc, char **argv) FROM " << __FILE__ << std::endl;
  init_gtest_tests(args);
  return run_gtest_tests();
}
""" % includes
    if includes:
        print('build_and_run_gtest_auto.py: making gtest_auto.gen.cpp')
        with open(cmake_dir + '/gtest_auto.gen.cpp', 'wb') as out:
            out.write(code)


if __name__ == '__main__':
    try:
        cmake_dir = get_cmake_dir('temp', cfg='Release')
    except AssertionError:
        rebuild_fast('gtest_wip')
        cmake_dir = get_cmake_dir('temp', cfg='Release')
    files = [x for x in sys.argv[1:] if x.endswith('.gtest.cpp')]
    files.extend(x.replace('.hpp', '.gtest.cpp') for x in sys.argv
                 if x.endswith('.hpp') and
                 os.path.exists(x.replace('.hpp', '.gtest.cpp')))
    make_gtest_auto_cpp(files, cmake_dir)
    assert not os.system('cd ' + cmake_dir + ' && ninja gtest_auto')
    assert not os.system(cmake_dir + '/gtest_auto')
