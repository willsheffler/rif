import shutil, os, glob, subprocess, re
from jinja2 import Template
from collections import namedtuple

# assume run from top project directory
pyfiles = subprocess.check_output('find src -regex [^.].+pybind.cpp'.split()).split()
pymodules = dict()
for f in pyfiles:
    if f == "src/riflib.pybind.cpp": continue
    print "sourcegen.py: found pybind file", f
    for line in subprocess.check_output(['grep','--with-filename','RIFLIB_PYBIND_', f]).splitlines():
        match = re.match("src/(.+).pybind.cpp:.* RIFLIB_PYBIND_(\w+)",line)
        path = match.group(1).replace('riflib/','')
        func = match.group(2)
        assert len(match.groups()) is 2
        pymodules[func] = path

# print pymodules

# todo: move the formatting into the jinja template, pass it the data not code!
forward = ''
code = ''

for v in set(pymodules.values()):
    code += '    py::module ' + v.replace('/','__') + ' = riflib.def_submodule("' + '").def_submodule("'.join(v.split('/')) + '");\n'
code += '\n'
for k,v in pymodules.items():
    forward += 'void RIFLIB_PYBIND_'+k+'(py::module & m);\n'
    code += '    RIFLIB_PYBIND_'+k+'('+v.replace('/','__')+');\n'

# print forward
# print code

if os.path.exists('src/riflib.pybind.cpp.TMP'):
    os.remove('src/riflib.pybind.cpp.TMP')

with open('src/riflib.pybind.cpp.jinja','r') as f:
    t = Template(f.read())
    with open('src/riflib.pybind.cpp.TMP','w') as o:
        newcontent = t.render(forward=forward, code=code)
        o.write(newcontent)

    diff = True
    if os.path.exists('src/riflib.pybind.cpp'):
        assert os.path.exists('src/riflib.pybind.cpp.TMP')
        diff = subprocess.call("diff src/riflib.pybind.cpp src/riflib.pybind.cpp.TMP".split())
    if diff:
        print 'sourcegen.py: updating riflib.pybind.cpp'
        if os.path.exists('src/riflib.pybind.cpp'):
            os.remove('src/riflib.pybind.cpp')
        os.rename('src/riflib.pybind.cpp.TMP','src/riflib.pybind.cpp')
    else:
        print "sourcegen.py: riflib.pybind.cpp is up to date"
        os.remove('src/riflib.pybind.cpp.TMP')
