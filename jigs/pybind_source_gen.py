#!/usr/bin/env python

"""generate main rifgen.gen.cpp pybind file with modules organized based on file paths
of *.pybind.cpp components"""

import subprocess
import os
import re
from collections import OrderedDict

import jinja2

def get_pybind_modules(srcpath):
    "find RIFLIB_PYBIND_ functions in *.pybind.cpp files"
    pbfiles = subprocess.check_output('find {} -regex [^.].+pybind.cpp'.format(srcpath).split())
    pymodules = OrderedDict()
    for pybindfile in pbfiles.splitlines():
        print "found pybind file", pybindfile
        try:
            grepped = subprocess.check_output(['grep', '-H', 'RIFLIB_PYBIND_', pybindfile])
        except:
            continue
        for line in grepped.splitlines():
            match = re.match(r".*src/(.+).pybind.cpp\:.* RIFLIB_PYBIND_(\w+)", line)
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
    diff = 1
    if os.path.exists(destfile):
        assert os.path.exists(testfile)
        diff = subprocess.call(['diff', testfile, destfile])
    if diff:
        print 'pybind_source_gen.py: updating', destfile
        if os.path.exists(destfile):
            os.remove(destfile)
        os.rename(testfile, destfile)
    else:
        print "pybind_source_gen.py: riflib.pybind.cpp is up to date"
        os.remove(testfile)

def shitty_make_code(pymodules):
    "todo: move the formatting into the jinja template, pass it the data not code!"
    code1 = ''
    code2 = ''
    for path in set(pymodules.values()):
        code2 += ('    py::module ' + path.replace('/', '__') +
                  ' = riflib.def_submodule("' +
                  '").def_submodule("'.join(path.split('/')) +
                  '");\n')
    code2 += '\n'
    for func, path in pymodules.items():
        code1 += 'void RIFLIB_PYBIND_' + func + '(py::module & m);\n'
        code2 += '    RIFLIB_PYBIND_' + func + '('+path.replace('/', '__') + ');\n'
    return code1, code2

def mkdir_if_necessary(path):
    os.system('mkdir -p '+path)

def mkfile_if_necessary(path, content):
    if not os.path.exists(path):
        print 'pybind_source_gen.py MAKING:', path
        with open(path,'w') as out:
            out.write(content)

def make_py_stencils(pymodules, srcdir):
    mkdir_if_necessary(srcdir)
    mkfile_if_necessary(srcdir+"/__init__.py",'from riflib_cpp import *\n'+
                        'import riflib_cpp\n'+
                        '__version__ = riflib_cpp.__version__\n')
    directories = set()
    for path in set(pymodules.values()):
        for nprefix in range(1,len(path.split('/'))):
            d = '/'.join(path.split('/')[:nprefix])
            directories.add(d)
            mkdir_if_necessary(srcdir + d)
    for path in directories:
        pyfile = srcdir + path + '/__init__.py'
        mkfile_if_necessary(pyfile, 'from riflib_cpp.'+path.replace('/','.')+' import *\n')
    for path in set(pymodules.values()):
        pyfile = srcdir+path+".py"
        mkfile_if_necessary(pyfile,'from riflib_cpp.'+path.replace('/','.')+' import *\n')

def main(template_fname, srcdir):
    "generate pybind sources"
    destfile = template_fname.replace('.jinja', '')
    pymodules = get_pybind_modules(srcdir) # assume in top level project dir
    make_py_stencils(pymodules, srcdir)
    forward, code = shitty_make_code(pymodules)
    with open(template_fname, 'r') as template_file:
        template = jinja2.Template(template_file.read())
    newcontent = template.render(forward=forward, code=code)
    update_file_if_needed(destfile, newcontent)

if __name__ == '__main__':
    import sys
    assert len(sys.argv) == 2
    srcdir = sys.argv[1]
    srcdir += '/' if not srcdir.endswith('/') else ''
    main(srcdir+'riflib.gen.cpp.jinja', srcdir)
