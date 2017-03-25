from __future__ import absolute_import, division, print_function
from builtins import bytes

import glob
import os
import sys
import multiprocessing
import re
import filecmp

import pytest

# todo this is a dup from setup.py....


def get_my_compiler():
    my_compiler = os.getenv('CXX', '').replace('/', '')
    if not my_compiler:
        my_compiler = "DEFAULT_CXX"
    return my_compiler


def get_my_python():
    return sys.executable.replace('/', '')


# todo: remove the above dups from setup.py

def in_ci_environment():
    return 'CI' in os.environ or 'READTHEDOCS' in os.environ


def get_proj_root():
    """find rood dir of this project by looking for .git and external"""
    proj_root = os.path.abspath('.')
    while proj_root != '/':
        if (not os.path.exists(proj_root + '/.git') and
                not os.path.exists(proj_root + '/external')):
            proj_root = '/'.join(proj_root.split('/')[:-1])
        else:
            break
    else:
        raise SystemError("can't find project root!")
    return proj_root


def get_build_dir(cfg='Release'):
    path = get_proj_root() + '/build_setup_py_' + cfg
    return path


def get_cmake_dir(prefix, cfg):
    """get directory setup.py builds stuff in, prefix is lib or temp"""
    version = '{}.{}'.format(sys.version_info.major, sys.version_info.minor)
    path = get_build_dir(cfg) + '/' + prefix + '*' + \
        version + '-' + get_my_python() + '-' + get_my_compiler()
    libdir = (glob.glob(path))
    assert len(libdir) == 1
    assert os.path.exists(libdir[0])
    return libdir[0]


def get_ignored_dirs(cfg):
    d, tgt = os.path.split(get_cmake_dir('temp', cfg))
    path = get_build_dir(cfg)
    decoys = os.listdir(path)
    # for x in decoys:
    # print(x)
    # print('tgt', tgt)
    assert tgt in decoys
    # print(len(decoys))
    decoys.remove(tgt)
    # print(len(decoys))
    assert not tgt in decoys
    return [d + '/' + x for x in decoys]


def which(program):
    """like linux which"""
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


def add_to_pypath(newpath):
    if isinstance(newpath, str):
        newpath = [newpath]
    if not hasattr(newpath, '__iter__'):
        newpath = [newpath]
    current = None
    if 'PYTHONPATH' in os.environ:
        current = os.environ['PYTHONPATH']
    os.environ['PYTHONPATH'] = ':'.join(os.path.abspath(p) for p in newpath)
    if current:
        os.environ['PYTHONPATH'] += ':' + current


def rebuild_setup_py_rif(cfg='Release'):
    proj_root = get_proj_root()
    assert not os.system('cd ' + proj_root + '; ' + sys.executable +
                         ' setup.py build --build-base=build_setup_py_' + cfg)


def make_docs(kind='html', cfg='Release'):
    proj_root = get_proj_root()
    rebuild_setup_py_rif()
    add_to_pypath(get_cmake_dir('lib', cfg=cfg))
    return os.system('cd ' + proj_root + '/docs; make ' + kind)


def rif_is_installed():
    rif_loader = False
    try:
        import pkgutil
        rif_loader = pkgutil.find_loader('rif')
    except ImportError:
        import importlib
        rif_loader = importlib.util.find_spec('rif')
    return rif_loader is not None


def is_devel_install():
    proj_root = get_proj_root()
    return proj_root + '/src' in sys.path


def remove_installed_rif():
    if rif_is_installed():
        os.system('echo y | ' + sys.executable + ' -m pip uninstall rif')
        assert is_devel_install() or not rif_is_installed()
        print('uninstalled rif from ' + sys.executable)


def get_ncpu():
    ncpu = multiprocessing.cpu_count()
    if ncpu > 4:
        ncpu = int(ncpu / 2)
    if in_ci_environment():
        ncpu = min(ncpu, 4)
        os.system('uname -a')
        print('build_utils.py: multiprocessing.cpu_count() = ',
              multiprocessing.cpu_count())
    return ncpu


def get_gtests(args_in):
    args = set(args_in)
    args.update(x.replace('.hpp', '.gtest.cpp') for x in args_in
                if os.path.exists(x.replace('.hpp', '.gtest.cpp')))
    gtests = set()
    for gtestfile in (x for x in args if x.endswith('.gtest.cpp')):
        print("    get_gtests", gtestfile)
        with open(gtestfile) as file:
            contents = file.read()
            for match in re.findall("TEST\(\s*(\S+?),\s*(\S+?)\s*\)", contents):
                assert len(match) is 2
                gtests.add(match[0])
    return gtests


def src_dir_new_file():
    fname = '.__build_utils_src_dir_contents'
    ftemp = '.__build_utils_src_dir_contents_tmp'
    os.system('find ./src -name \\*.cpp > ' + ftemp)
    are_same = filecmp.cmp(fname, ftemp) if os.path.exists(fname) else False
    if not are_same:
        os.rename(ftemp, fname)
    else:
        os.remove(ftemp)
    return not are_same


def rebuild_fast(target='rif_cpp', cfg='Release', force_redo_cmake=False):
    try:
        if force_redo_cmake:
            raise Exception
        cmake_dir = get_cmake_dir('temp', cfg=cfg)
        makeexe = 'ninja'
        if not which('ninja'):
            makeexe = 'make'
        ncpu = multiprocessing.cpu_count()
        return os.system('cd ' + cmake_dir + '; ' + makeexe + ' -j%i ' % ncpu + target)
    except Exception as e:
        return rebuild_setup_py_rif(cfg=cfg)


def build_and_test():
    assert os.path.exists('CMakeLists.txt') and os.path.exists('setup.py')
    print("== build_and_test ==")
    cfg = 'Release'
    testfiles = [x for x in sys.argv[1:] if x.endswith('.py') and
                 os.path.basename(x).startswith('test')]
    pybindfiles = [x for x in sys.argv[1:] if x.endswith('.pybind.cpp')]
    gtests = get_gtests(sys.argv[1:])
    print('calling rebuild_fast')
    no_xdist = len(testfiles) or len(gtests)
    no_xdist |= os.system('egrep "#.*cpp_files" pytest.ini') == 0
    force_redo_cmake = len(pybindfiles) or not no_xdist or src_dir_new_file()
    rebuild_fast(target='rif_cpp gtest_all',
                 cfg=cfg, force_redo_cmake=force_redo_cmake)

    # TODO both here and in docs, this gets messed
    #      up when rif is actually installed
    libdir = os.path.abspath(get_cmake_dir('lib', cfg))
    print('adding to python path:', libdir)
    assert os.path.exists(libdir)
    # need to use sys.path for this process
    sys.path.append(libdir)
    # need to use PYTHONPATH env for xdist subprocesses
    add_to_pypath(libdir)
    assert libdir in sys.path
    assert libdir in os.environ['PYTHONPATH'].split(':')

    ncpu = get_ncpu()
    args = list(testfiles)
    if not (testfiles or gtests):
        args = ['.']
    if not no_xdist:
        args.extend('--cov=./src -n{}'.format(ncpu).split())
    if gtests:
        args.extend(['-k', ' or '.join(gtests)])
    args.extend('--ignore build'.split())
    if testfiles and not gtests:
        args.extend('--ignore build_setup_py_Release'.split())
    for decoy in get_ignored_dirs(cfg):
        args += ['--ignore', decoy]
    print('==================================================================================')
    sys.stdout.write('pytest ')
    for arg in args:
        if arg.startswith('-'):
            sys.stdout.write('\n   ')
        sys.stdout.write(' ' + arg)
    print()
    sys.argv[1:] = args
    assert not pytest.main()




def make_gtest_auto_cpp(files, cmake_dir):
    includes = "\n".join("#include <%s>" % f for f in files)
    code = """#include "test/gtest_util.hpp"
%s
int main(int argc, char **argv) {
  std::vector<std::string> args;
  for (int i = 0; i < argc; ++i) args.push_back(std::string(argv[i]));
  std::cout << "int main(int argc, char **argv) FROM " << __FILE__ << std::endl;
  init_gtest_tests(args);
  return run_gtest_tests();
}
""" % includes
    if includes:
        print('build_and_run_gtest_auto.py: making gtest_auto.gen.cpp')
        with open(cmake_dir + '/gtest_auto.gen.cpp', 'wb') as out:
            out.write(code)



def build_and_run_gtest_auto():

    try:
        cmake_dir = get_cmake_dir('temp', cfg='Release')
    except AssertionError:
        rebuild_fast('gtest_wip')
        cmake_dir = get_cmake_dir('temp', cfg='Release')
    files = [x for x in sys.argv[1:] if x.endswith('.gtest.cpp')]
    files.extend(x.replace('.hpp', '.gtest.cpp') for x in sys.argv
                 if x.endswith('.hpp') and
                 os.path.exists(x.replace('.hpp', '.gtest.cpp')))
    if len(files) + 1 is not len(sys.argv):
        raise NotImplementedError
    make_gtest_auto_cpp(files, cmake_dir)
    assert not os.system('cd ' + cmake_dir + ' && ninja gtest_auto')
    assert not os.system(cmake_dir + '/gtest_auto')
