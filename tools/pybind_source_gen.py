#!/usr/bin/env python

from __future__ import print_function

"""generate main rifgen.gen.cpp pybind file with modules
organized based on file paths of *.pybind.cpp components"""

import subprocess
import os
import re
import datetime
from collections import OrderedDict
import sys

try:
    import jinja2

except ImportError as error:
    print('!!!!!!!!!!!!!! ' * 100)
    print('cant import jinja2')
    print(sys.executable)
    print('SYS.PATH')
    for p in sys.path:
        print('   ', p)
        try:
            for d in os.listdir(p):
                if 'jinja' in d:
                    print('       ', d)
        except OSError:
            print('cant read', d)
    sys.exit(-1)


def get_pybind_modules(srcpath):
    "find RIFLIB_PYBIND_ functions in *.pybind.cpp files"
    pymodules = OrderedDict()
    for root, _, files in os.walk(srcpath):
        for basename in (x for x in files if x.endswith('.pybind.cpp')):
            pybindfile = root + '/' + basename
            print("found pybind file", pybindfile)
            try:
                # todo: replace this with python
                try:
                    grepped = subprocess.check_output(
                        ['grep', '-H', 'RIFLIB_PYBIND_', pybindfile])
                except subprocess.CalledProcessError:
                    grepped = ''
                # print('grepped', grepped)
            except OSError:
                continue
            for line in grepped.splitlines():
                line = str(line)
                match = re.match(
                    r".*src/rif/(.+).pybind.cpp\:.* RIFLIB_PYBIND_(\w+)", line)
                # print 'line:', line
                # print match
                # assert len(match.groups()) is 2
                path = match.group(1)
                func = match.group(2)
                pymodules[func] = path
    return pymodules


def update_file_if_needed(destfile, newcontent):
    "update destfile with newcontent if needed, otherwise don't modifiy file"
    testfile = destfile + ".TMP"
    if os.path.exists(testfile):
        os.remove(testfile)
    with open(testfile, 'w') as out:
        out.write(newcontent)
    diff = ''
    if os.path.exists(destfile):
        assert os.path.exists(testfile)
        try:
            need to get actual diff from this
            then figure out why pybind def_submodules are messed
            then why can't make a SceneBase
            diff = subprocess.call(['diff', testfile, destfile])
        except subprocess.CalledProcessError:
            pass
        print('diff:', diff)
        if len(diff.splitlines()) is 4 and 'compiled on' in diff:
            diff = None
    if diff:
        print('pybind_source_gen.py: updating', destfile)
        if os.path.exists(destfile):
            os.remove(destfile)
        os.rename(testfile, destfile)
    else:
        print("pybind_source_gen.py: rif.pybind.cpp is up to date")
        os.remove(testfile)


def shitty_make_code(pymodules):
    """todo: move the formatting into the jinja
    template, pass it the data not code!"""
    code1 = ''
    code2 = ''
    for path in set(pymodules.values()):
        code2 += ('    py::module ' + path.replace('/', '__') +
                  ' = rif.def_submodule("' +
                  '").def_submodule("'.join(path.split('/')) +
                  '");\n')
    code2 += '\n'
    for func, path in list(pymodules.items()):
        code1 += 'void RIFLIB_PYBIND_' + func + '(py::module & m);\n'
        code2 += '    RIFLIB_PYBIND_' + func + \
            '(' + path.replace('/', '__') + ');\n'
    return code1, code2


def mkdir_if_necessary(path):
    os.system('mkdir -p ' + path)


def mkfile_if_necessary(path, content):
    if not os.path.exists(path):
        print('pybind_source_gen.py MAKING:', path)
        with open(path, 'w') as out:
            out.write(content)


def make_py_stencils(pymodules, srcdir):
    mkdir_if_necessary(srcdir)
    mkfile_if_necessary(srcdir + "/__init__.py", 'from rif_cpp import *\n' +
                        'import rif_cpp\n' +
                        '__version__ = rif_cpp.__version__\n')
    directories = set()
    for path in set(pymodules.values()):
        for nprefix in range(1, len(path.split('/'))):
            d = '/'.join(path.split('/')[:nprefix])
            directories.add(d)
            mkdir_if_necessary(srcdir + d)
    for path in directories:
        pyfile = srcdir + path + '/__init__.py'
        mkfile_if_necessary(pyfile, 'from rif_cpp.' +
                            path.replace('/', '.') + ' import *\n')
    for path in set(pymodules.values()):
        pyfile = srcdir + path + ".py"
        mkfile_if_necessary(pyfile, 'from rif_cpp.' +
                            path.replace('/', '.') + ' import *\n')


def main(template_fname, srcdir, dstdir):
    "generate pybind sources"
    destfile = dstdir + os.path.basename(template_fname.replace('.jinja', ''))
    pymodules = get_pybind_modules(srcdir)  # assume in top level project dir
    make_py_stencils(pymodules, srcdir)
    forward, code = shitty_make_code(pymodules)
    with open(template_fname, 'r') as template_file:
        template = jinja2.Template(template_file.read())
    newcontent = template.render(
        forward=forward, code=code, mydate=str(datetime.datetime.now()))
    update_file_if_needed(destfile, newcontent)

if __name__ == '__main__':
    print("== pybind_source_gen.py ==")
    import sys
    assert len(sys.argv) == 3
    srcdir = sys.argv[1]
    dstdir = sys.argv[2]
    srcdir += '/' if not srcdir.endswith('/') else ''
    dstdir += '/' if not dstdir.endswith('/') else ''
    main(srcdir + 'rif.gen.cpp.jinja', srcdir, dstdir)
