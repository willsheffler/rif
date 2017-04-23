#!/usr/bin/env python
from __future__ import print_function

import os
import sys
import subprocess
import shutil


def main():
    src = sys.argv[1]
    dst = sys.argv[2]
    exts = sys.argv[3:]
    assert exts
    assert src != dst
    for root, dirs, files in os.walk(src):
        for file in files:
            if any(file.endswith(ext) for ext in exts):
                newfile = root.replace(src, dst) + '/' + file
                os.system('mkdir -p ' + os.path.dirname(newfile))
                shutil.copy(root + '/' + file, newfile)


if __name__ == '__main__':
    main()
