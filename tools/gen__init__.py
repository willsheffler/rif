import sys
import os
import subprocess


def all_parent_dirs(path):
    spath = path.split('/')
    return ['/'.join(spath[:i]) for i in range(1, len(spath))]


def main():
    assert len(sys.argv) is 2
    src = sys.argv[1]
    pyfiles = subprocess.check_output(
        'find {} -regex [^.].+[.]py'.format(src).split()).split()
    packages = set()
    for pyfile in pyfiles:
        if os.path.basename(pyfile).startswith('_test'):
            continue
        packages.update(all_parent_dirs(pyfile))
    for p in packages:
        if not os.path.exists(p + '/__init__.py'):
            print 'creating empty', p + '/__init__.py'
            with open(p + '/__init__.py', 'w'):
                pass


if __name__ == '__main__':
    main()
