import sys
# import os
import subprocess


def main():
    assert len(sys.argv) is 2
    src = sys.argv[1]
    pyfiles = subprocess.check_output(
        'find {} -regex [^.].+[.]py'.format(src).split())
    for pyfile in pyfiles:
        print 'PYFILE', pyfile

if __name__ == '__main__':
    main()
