"""generate empty __init__.py packages files where necessary"""

from __future__ import print_function

import sys
import os


def all_parent_dirs(path, prefix):
    spath = path.split('/')
    # todo have it check against the arg instead of -1 below
    allparents = ['/'.join(spath[:i]) for i in range(1, len(spath))]
    return [x for x in allparents if x.startswith(prefix)]


def main():
    print("== gen__init__.py starting ==")
    assert len(sys.argv) is 2
    src = sys.argv[1]
    packages = set()
    extensions = '.py .pybind.cpp'.split()
    for root, _, files in os.walk(src):
        for fname in files:
            if (any(fname.endswith(x) for x in extensions) and
                not (fname.startswith('test') or
                     fname.endswith('test.py'))):
                packages.update(all_parent_dirs(root + '/dummy', src))
    for p in packages:
        if not os.path.exists(p + '/__init__.py'):
            print('creating empty', p + '/__init__.py')
            with open(p + '/__init__.py', 'w') as out:
                out.write('"""docstring for ' + p.replace('src/rif/', 'rif/') +
                          '/__init__.py\n"""')


if __name__ == '__main__':
    main()
