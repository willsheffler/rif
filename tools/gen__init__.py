import sys
import os
import subprocess


def all_parent_dirs(path):
    spath = path.split('/')
    return ['/'.join(spath[:i]) for i in range(1, len(spath))]


def main():
    assert len(sys.argv) is 2
    src = sys.argv[1]
    # pyfiles = subprocess.check_output(
        # 'find {} -regex [^.].+[.]py'.format(src).split()).split()
    packages = set()
    for root, _, files in os.walk(src):


        for fname in files:
            if (not fname.endswith('.py') or fname.startswith('test') or 
                fname.endswith('test.py')):
                continue
            packages.update(all_parent_dirs(root+'/dummy'))
    for p in packages:
        if not os.path.exists(p + '/__init__.py'):
            print('creating empty', p + '/__init__.py')
            with open(p + '/__init__.py', 'w'):
                pass


if __name__ == '__main__':
    main()
