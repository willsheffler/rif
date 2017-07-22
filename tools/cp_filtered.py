#!/usr/bin/env python
import os
import sys
import shutil
try:
    import yaml
    _have_yaml = True
except ImportError:
    _have_yaml = False


def main():
    src = sys.argv[1]
    dst = sys.argv[2]
    exts = sys.argv[3:]
    assert exts
    assert src != dst
    if _have_yaml:
        exclude = set(yaml.load(os.path.join(src, 'do_not_package.yml')))
    else:
        exclude = []
    for root, dirs, files in os.walk(src):
        for file in files:
            if file in exclude or os.path.basename(file) in exclude:
                continue
            if any(file.endswith(ext) for ext in exts):
                newfile = root.replace(src, dst) + '/' + file
                os.system('mkdir -p ' + os.path.dirname(newfile))
                shutil.copy(root + '/' + file, newfile)


if __name__ == '__main__':
    main()
