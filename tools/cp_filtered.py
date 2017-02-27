#!/usr/bin/env python
from __future__ import print_function

import os
import sys
import subprocess
import shutil


def main():
    ext = sys.argv[1]
    src = sys.argv[2]
    dst = sys.argv[3]
    assert ext
    assert src != dst
    for root, dirs, files in os.walk(src):
        for file in files:
            if file.endswith(ext):
                newfile = root.replace(src, dst) + '/' + file
                os.system('mkdir -p ' + os.path.dirname(newfile))
                shutil.copy(root + '/' + file, newfile)


if __name__ == '__main__':
    main()
