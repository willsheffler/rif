from __future__ import absolute_import, division, print_function

import os
import re
import sys
import platform
import subprocess
import multiprocessing
from glob import glob
from collections import defaultdict

from setuptools import setup, Extension, find_packages
from setuptools.command.build_ext import build_ext
from distutils.version import LooseVersion

# todo: have these utility functions in a module... but how to import them
# into setup.py...


def get_my_compiler():
    my_compiler = os.getenv('CXX', '').replace('/', '')
    if not my_compiler:
        my_compiler = "DEFAULT_CXX"
    return my_compiler


def get_my_python():
    return sys.executable.replace('/', '')


def my_getenv(name):
    if name in os.environ:
        return os.environ[name]
    else:
        return "DEFAULT_" + name


def in_conda():
    return ('Anaconda' in sys.version or
            'Continuum Analytics' in sys.version or
            'conda' in sys.executable
            )


def which(program):
    import os

    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            path = path.strip('"')
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file

    return None


def infer_config_from_build_dirname(path):
    path = path.split('/')[0]
    if path.startswith('build_setup_py_'):
        return path.replace('build_setup_py_', '')


_rif_setup_opts = defaultdict(list)
_remove_from_sys_argv = list()
for arg in sys.argv:
    if arg.startswith('--rif_setup_opts_'):
        flag, val = arg.split('=')
        _rif_setup_opts[flag[17:]] = val.split(',')
        _remove_from_sys_argv.append(arg)
for arg in _remove_from_sys_argv:
    sys.argv.remove(arg)
print('setup.py: rif args:')
for flag, val in _rif_setup_opts.items():
    print('    ', flag, '=', val)

print('setup.py: compiler:', get_my_compiler())
print('setup.py: python:', get_my_python())
for evar in "CC CXX CXXFLAGS".split():
    print('sepup.py: env:', evar, my_getenv(evar))


class CMakeExtension(Extension):

    def __init__(self, name, sourcedir=''):
        Extension.__init__(self, name, sources=[])
        self.sourcedir = os.path.abspath(sourcedir)


class CMakeBuild(build_ext):

    def __init__(self, *args, **kwargs):
        build_ext.__init__(self, *args, **kwargs)
        self.my_tag = get_my_python() + '-' + get_my_compiler()

    def run(self):
        if not self.build_temp.endswith(self.my_tag):
            self.build_temp += '-' + self.my_tag
        try:
            out = subprocess.check_output(['cmake', '--version'])
        except OSError:
            raise RuntimeError(
                "CMake must be installed to build the following extensions: " +
                ", ".join(e.name for e in self.extensions))

        if platform.system() == "Windows":
            cmake_version = LooseVersion(
                re.search(r'version\s*([\d.]+)', out.decode()).group(1))
            if cmake_version < '3.1.0':
                raise RuntimeError("CMake >= 3.1.0 is required on Windows")

        for ext in self.extensions:
            self.build_extension(ext)

    def get_ext_fullpath(self, ext_name):
        defaultname = build_ext.get_ext_fullpath(self, ext_name)
        path = os.path.dirname(defaultname)
        path += '-' + self.my_tag + '/'
        path += os.path.basename(defaultname)
        return path, defaultname

    def build_extension(self, ext):
        extdir, defaultextdir = self.get_ext_fullpath(ext.name)
        extdir = os.path.abspath(os.path.dirname(extdir))
        defaultextdir = os.path.abspath(os.path.dirname(defaultextdir))
        cmake_args = ['-DCMAKE_LIBRARY_OUTPUT_DIRECTORY=' + extdir,
                      '-DPYTHON_EXECUTABLE=' + sys.executable,
                      ]

        # if 'CI' in os.environ:
        #     print('checking for CMAKE_OPTIONS')
        #     for k, v in os.environ.items():
        #         print('ENV', k, '=', v)
        #         if k is 'CMAKE_OPTIONS':
        #             print('finding any *.so in local boost root')
        #             os.system('find %s -name \*.so' %
        #                       v.replace('-DBOOST_ROOT=', ''))
        if 'CMAKE_OPTIONS' in os.environ:
            print('setup.py add CMAKE_OPTIONS:', os.environ['CMAKE_OPTIONS'])
            cmake_args += os.environ['CMAKE_OPTIONS'].split()
        ncpu = multiprocessing.cpu_count()
        if 'CI' in os.environ or 'READTHEDOCS' in os.environ:
            ncpu = 4
            os.system('uname -a')
            print('setup.py: multiprocessing.cpu_count() is ',
                  multiprocessing.cpu_count())
        if which('ninja'):
            cmake_args.append('-GNinja')
        cfg = infer_config_from_build_dirname(self.build_temp)
        if not cfg:
            cfg = 'Debug' if self.debug else 'Release'
        build_args = ['--config', cfg]

        if platform.system() == "Windows":
            cmake_args += [
                '-DCMAKE_LIBRARY_OUTPUT_DIRECTORY_{}={}'.format(cfg.upper(), extdir)]
            if sys.maxsize > 2**32:
                cmake_args += ['-A', 'x64']
                build_args += ['--', '/m']
        else:
            cmake_args += ['-DCMAKE_BUILD_TYPE=' + cfg]
            build_args += ['--', '-j' + str(ncpu)]
        cmake_args += _rif_setup_opts['cmake_args']
        build_args += _rif_setup_opts['build_args']
        cmake_args += ['-DPYTHON_EXECUTABLE:FILEPATH=' + sys.executable, ]

        env = os.environ.copy()
        env['CXXFLAGS'] = '{} -DVERSION_INFO=\\"{}\\"'.format(
            env.get('CXXFLAGS', ''), self.distribution.get_version())
        if in_conda():
            condadir = os.path.dirname(sys.executable)[:-4]
            env['CXXFLAGS'] = env['CXXFLAGS'] + ' -I' + condadir + '/include'
            env['CXXFLAGS'] = env['CXXFLAGS'] + ' -L' + condadir + '/lib'
            print('setup.py: adding -I/-L for conda', condadir)
        if not os.path.exists(self.build_temp):
            os.makedirs(self.build_temp)
        try:
            subprocess.check_call(
                ['cmake', ext.sourcedir] + cmake_args, cwd=self.build_temp, env=env)
            subprocess.check_call(
                ['cmake', '--build', '.'] + build_args, cwd=self.build_temp)
        except subprocess.CalledProcessError as e:
            print('setup.py exiting with returncode', e.returncode)
            sys.exit(e.returncode)
            print('this should never happen')
        if os.path.exists(defaultextdir):
            os.remove(defaultextdir)
        # TODO: figure out how to remove this, needed for tox
        os.symlink(extdir, defaultextdir)


setup(
    name='rif',
    version='0.0.1',
    author='Will Sheffler',
    author_email='willsheffler@gmail.com',
    description='Rotamer Interaction Field protein design library',
    long_description='',
    url='https://github.com/willsheffler/rif',
    ext_modules=[CMakeExtension('rif_cpp')],
    cmdclass=dict(build_ext=CMakeBuild),
    zip_safe=False,
    setup_requires=['wheel', 'pytest-runner'],
    tests_require=['pytest', 'pytest-xdist', 'tox',
                   'hypothesis', 'colorama', 'pytest_cpp', 'jinja2'],
    test_suite='pytest',
)
