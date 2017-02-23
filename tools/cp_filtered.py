#!/usr/bin/env python

import os
import sys
import subprocess
import shutil


def main():
    ext = sys.argv[1]
    dst = sys.argv[3]
    src = sys.argv[2]
    extfiles = subprocess.check_output(
        'find {} -regex [^.].+{}'.format(src, ext).split())
    assert extfiles
    for extfile in extfiles.splitlines():
        newfile = extfile.replace(src.encode('utf-8'), dst.encode('utf-8'))
        os.system(b'mkdir -p ' + os.path.dirname(newfile))
        shutil.copy(extfile, newfile)


if __name__ == '__main__':
    main()
