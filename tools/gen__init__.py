import sys
import os
import subprocess


def all_parent_dirs(path):
    spath = path.split(b'/')
    return [b'/'.join(spath[:i]) for i in range(1, len(spath))]


def main():
    assert len(sys.argv) is 2
    src = sys.argv[1]
    pyfiles = subprocess.check_output(
        'find {} -regex [^.].+[.]py'.format(src).split()).split()
    packages = set()
    for pyfile in pyfiles:
        fname = str(os.path.basename(pyfile))
        if (fname.startswith('test') or fname.endswith('test')):
            continue
        packages.update(all_parent_dirs(pyfile))
    for p in packages:
        if not os.path.exists(p + b'/__init__.py'):
            print('creating empty', p + b'/__init__.py')
            with open(p + b'/__init__.py', 'w'):
                pass


if __name__ == '__main__':
    main()
