#!/usr/bin/env python
import subprocess
import sys
import os


def main():
    # print os.getcwd()
    if len(sys.argv) is 1:
        wd = '.'
    elif len(sys.argv) is 2:
        wd = sys.argv[1]
    else:
        raise ValueError
    count = 0
    for root, dirs, files in os.walk(wd):
        for file in files:
            path = os.path.join(root, file)
            if file.endswith('.pyc'):
                # print 'removing', path
                os.remove(path)
                count += 1
            if 'build_setup_py' in root and file.endswith('.py'):
                # print 'removing', path
                os.remove(path)
                count += 1
    print('removed', count, '.pyc and build/*.py files')


if __name__ == '__main__':
    main()
