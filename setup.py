import os
import re
import sys
import platform
import subprocess
import multiprocessing
from glob import glob
from collections import defaultdict

# assert sys.version_info.major > 2

# import numpy

from setuptools import setup, Extension, find_packages
from setuptools.command.build_ext import build_ext
from distutils.version import LooseVersion


setup_requires = ['wheel', 'pytest-runner', 'jinja2']
install_requires = ['numpy', 'pandas', 'xarray', 'parsimonious', 'mock', 'pytest',
                    'pytest-xdist', 'hypothesis', 'colorama', 'pytest_cpp',
                    'pytest-sugar', 'tqdm', 'dask', 'homog', 'jinja2']
tests_require = ['pytest', 'pytest-xdist', 'hypothesis', 'colorama',
                 'pytest_cpp', 'homog', 'jinja2']


###############################################################################
# horrible duplicates with setup.py and tools/build_utils.py
# setup.py must stand alone for detox to work(?)
# and importing from setup.py is problematic....
###############################################################################

def get_my_compiler():
    my_compiler = os.getenv('CXX', '').replace('/', '')
    if not my_compiler:
        my_compiler = "CXX"
    my_compiler = my_compiler.replace('usrbin', 'UB')
    my_compiler = my_compiler.replace('clang++', 'C')
    my_compiler = my_compiler.replace('g++', 'G')
    my_compiler = my_compiler.replace('-', '')
    my_compiler = my_compiler.replace('.', '')
    return my_compiler


def get_my_python():
    home = os.environ['HOME'] if 'HOME' in os.environ else "HOME"
    my_python = sys.executable
    my_python = my_python.replace('/usr/bin', 'UB')
    my_python = my_python.replace('/bin', 'B')
    my_python = my_python.replace('.tox', 'T')
    my_python = my_python.replace('anaconda', 'A')
    my_python = my_python.replace('miniconda', 'A')
    my_python = my_python.replace('conda/envs', 'ce')
    my_python = my_python.replace('software', 'S')
    my_python = my_python.replace(home, 'H')
    my_python = my_python.replace('python', 'Py')
    my_python = my_python.replace('rif', 'R')
    my_python = my_python.replace('/', '')
    return my_python


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
    path = os.path.basename(os.path.dirname(path))
    if path is 'buildD':
        return "Debug"
    return 'Release'

###############################################################################


_rif_setup_opts = defaultdict(list)
_remove_from_sys_argv = list()
for arg in sys.argv:
    if arg.startswith('--rif_setup_opts_'):
        flag, val = arg.split('=')
        _rif_setup_opts[flag[17:]] = val.split(',')
        _remove_from_sys_argv.append(arg)
for arg in _remove_from_sys_argv:
    sys.argv.remove(arg)
if _rif_setup_opts:
    print('setup.py: rif args:')
    for flag, val in _rif_setup_opts.items():
        print('    ', flag, '=', val)

# print('setup.py: compiler:', get_my_compiler())
# print('setup.py: python:', get_my_python())
# for evar in "CC CXX CXXFLAGS".split():
    # print('setup.py: env:', evar, my_getenv(evar))


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
            if 'READTHEDOCS' in os.environ:
                ncpu = 1
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
                '-DCMAKE_LIBRARY_OUTPUT_DIRECTORY_{}={}'.format(cfg.upper(),
                                                                extdir)]
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
        # env['CXXFLAGS'] += ' -I' + numpy.get_include()
        if in_conda():
            condadir = os.path.dirname(sys.executable)[:-4]  # /bin
            env['CXXFLAGS'] = env['CXXFLAGS'] + ' -I' + condadir + '/include'
            # this causes build failure....
            # env['CXXFLAGS'] = env['CXXFLAGS'] + ' -L' + condadir + '/lib'
            print('setup.py: adding -I/-L for conda', condadir)
        if not os.path.exists(self.build_temp):
            os.makedirs(self.build_temp)
        try:
            try:
                with open('log/cmake_cmd.log', 'w') as out:
                    out.writelines([
                        'setup.py: ' + os.linesep,
                        'CXXFLAGS: ' + env['CXXFLAGS'],
                        'cmake ' + ext.sourcedir + ' ' + ' '.join(cmake_args),
                        'cmake --build . ' + ' '.join(build_args),
                    ])
            except Exception:
                pass
            subprocess.check_call(
                ['cmake', ext.sourcedir] + cmake_args, cwd=self.build_temp, env=env)
            print("setup.py: cmake --build", " ".join(cmake_args))
            subprocess.check_call(
                ['cmake', '--build', '.'] + build_args, cwd=self.build_temp)
        except subprocess.CalledProcessError as e:
            print('setup.py exiting with returncode', e.returncode)
            sys.exit(e.returncode)
            print('this should never happen')
        # TODO: figure out how to remove this, needed for tox
        if os.path.exists(defaultextdir):
            os.remove(defaultextdir)
        os.symlink(extdir, defaultextdir)


# took these out... goal: have cmake manage everything
# def isdatfile(f):
#     return f.endswith('.gz') or f.endswith('csv') or f.endswith('.dat')


# def marshal_package_data():
#     pdat = list()
#     for root, dirs, files in os.walk('src/rif'):
#         d = root.replace('src/rif/', '')
#         datfiles = [os.path.join(d, f) for f in files if isdatfile(f)]
#         if datfiles:
#             print('marshaling datfiles:', d)
#             pdat.extend(datfiles)
#     return pdat


setup(
    name='rif',
    version='0.0.1',
    author='Will Sheffler',
    author_email='willsheffler@gmail.com',
    description='Rotamer Interaction Field protein design library',
    long_description='',
    url='https://github.com/willsheffler/rif',
    ext_modules=[CMakeExtension('_rif')],
    cmdclass=dict(build_ext=CMakeBuild),
    zip_safe=False,
    setup_requires=setup_requires,
    install_requires=install_requires,
    tests_require=tests_require,
    test_suite='pytest',
    # packages=['rif'],
    # package_dir={'rif': 'src/rif'},
    # package_data={'rif': marshal_package_data()}
)
