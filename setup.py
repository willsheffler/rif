from __future__ import absolute_import, division, print_function

import os
import re
import sys
import platform
import subprocess
import multiprocessing
from glob import glob

from setuptools import setup, Extension, find_packages
from setuptools.command.build_ext import build_ext
from distutils.version import LooseVersion


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


def infer_config(path):
    path = path.split('/')[0]
    if path.startswith('build_setup_py_'):
        return path.replace('build_setup_py_', '')


class CMakeExtension(Extension):

    def __init__(self, name, sourcedir=''):
        Extension.__init__(self, name, sources=[])
        self.sourcedir = os.path.abspath(sourcedir)


class CMakeBuild(build_ext):

    def run(self):
        try:
            out = subprocess.check_output(['cmake', '--version'])
        except OSError:
            raise RuntimeError("CMake must be installed to build the following extensions: " +
                               ", ".join(e.name for e in self.extensions))

        if platform.system() == "Windows":
            cmake_version = LooseVersion(
                re.search(r'version\s*([\d.]+)', out.decode()).group(1))
            if cmake_version < '3.1.0':
                raise RuntimeError("CMake >= 3.1.0 is required on Windows")

        for ext in self.extensions:
            self.build_extension(ext)

    def build_extension(self, ext):
        extdir = os.path.abspath(os.path.dirname(
            self.get_ext_fullpath(ext.name)))
        cmake_args = ['-DCMAKE_LIBRARY_OUTPUT_DIRECTORY=' + extdir,
                      '-DPYTHON_EXECUTABLE=' + sys.executable,
                      ]
        if which('ninja'):
            cmake_args.append('-GNinja')

        cfg = infer_config(self.build_temp)
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
            # build_args += ['--', '-j'+str(multiprocessing.cpu_count()), 'rif']
            build_args += ['--', '-j' + str(multiprocessing.cpu_count())]

        env = os.environ.copy()
        env['CXXFLAGS'] = '{} -DVERSION_INFO=\\"{}\\"'.format(env.get('CXXFLAGS', ''),
                                                              self.distribution.get_version())
        if not os.path.exists(self.build_temp):
            os.makedirs(self.build_temp)
        subprocess.check_call(['cmake', ext.sourcedir] +
                              cmake_args, cwd=self.build_temp, env=env)
        subprocess.check_call(['cmake', '--build', '.'] +
                              build_args, cwd=self.build_temp)


setup(
    name='rif',
    version='0.0.1',
    author='Will Sheffler',
    author_email='willsheffler@gmail.com',
    description='C++ based support libraries for rif',
    long_description='',
    ext_modules=[CMakeExtension('rif')],
    cmdclass=dict(build_ext=CMakeBuild),
    zip_safe=False,
    test_suite='pytest',
    tests_require=['pytest', 'hypothesis', 'colorama', 'pytest_cpp', 'jinja2'],
)
