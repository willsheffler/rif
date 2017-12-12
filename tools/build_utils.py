# from builtins import bytes

import glob
import os
from os.path import join, abspath, dirname, basename, exists
import sys
import multiprocessing
import shutil
import re
import filecmp
print(sys.version)
for p in sys.path:
    print(" ", p)
import pytest

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
            exe_file = join(path, program)
            if is_exe(exe_file):
                return exe_file
    return None


def infer_config_from_build_dirname(path):
    path = basename(dirname(path))
    if path is 'buildD':
        return "Debug"
    return 'Release'

###############################################################################


def get_build_dir(cfg='Release'):
    if cfg is 'Release':
        d = 'buildR'
    elif cfg is 'Debug':
        d = 'buildD'
    else:
        raise NotImplemented('cfg not recognised: ' + cfg)
    return join(get_proj_root(), d)


def in_CI_environment():
    return 'CI' in os.environ or 'READTHEDOCS' in os.environ


_VERBOSE = in_CI_environment()


def get_proj_root():
    """find rood dir of this project by looking for .git and external"""
    proj_root = abspath('.')
    while proj_root != '/':
        if (not exists(proj_root + '/.git') and
                not exists(proj_root + '/external')):
            proj_root = '/'.join(proj_root.split('/')[:-1])
        else:
            break
    else:
        raise SystemError("can't find project root!")
    return proj_root


def get_cmake_dir(prefix, cfg):
    """get directory setup.py builds stuff in, prefix is lib or temp"""
    version = '{}.{}'.format(sys.version_info.major, sys.version_info.minor)
    path = join(get_build_dir(cfg), prefix + '*' + version + '-' +
                get_my_python() + '-' + get_my_compiler())
    libdir = glob.glob(path)
    if not libdir:
        return None
    if len(libdir) is not 1:
        print('    libdir:', libdir)
    assert len(libdir) is 1
    assert exists(libdir[0])
    return libdir[0]


def get_ignored_dirs(cfg):
    d, lib = os.path.split(get_cmake_dir('lib', cfg))
    d, tgt = os.path.split(get_cmake_dir('temp', cfg))
    path = get_build_dir(cfg)
    decoys = os.listdir(path)
    # for x in decoys:
    # print(x)
    # print('tgt', tgt)
    assert tgt in decoys
    assert lib in decoys
    # print(len(decoys))
    decoys.remove(tgt)
    decoys.remove(lib)
    # print(len(decoys))
    assert tgt not in decoys
    return [join(d, x) for x in decoys]


def add_to_pypath(newpath):
    if isinstance(newpath, str):
        newpath = [newpath]
    if not hasattr(newpath, '__iter__'):
        newpath = [newpath]
    current = None
    if 'PYTHONPATH' in os.environ:
        current = os.environ['PYTHONPATH']
    os.environ['PYTHONPATH'] = ':'.join(abspath(p) for p in newpath)
    if _VERBOSE:
        print('== added to PYTHONPATH:', os.environ['PYTHONPATH'], '==')
    if current:
        os.environ['PYTHONPATH'] += ':' + current


def rebuild_rif(cfg='Release'):
    proj_root = get_proj_root()
    build_dir = get_build_dir(cfg)
    assert not os.system('cd ' + proj_root + '; ' + sys.executable +
                         ' setup.py build --build-base=' + build_dir)


def make_docs(kind='html', cfg='Release'):
    proj_root = get_proj_root()
    rebuild_rif()
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
        ncpu = int(ncpu * 4.0 / 8.0)
    if in_CI_environment():
        ncpu = min(ncpu, 4)
        os.system('uname -a')
        print('build_utils.py: multiprocessing.cpu_count() = ',
              multiprocessing.cpu_count())
    return ncpu


def get_gtests(args_in):
    args = set(args_in)
    args.update(x.replace('.hpp', '.gtest.cpp') for x in args_in
                if exists(x.replace('.hpp', '.gtest.cpp')))
    gtests = set()
    for gtestfile in (x for x in args if x.endswith('.gtest.cpp')):
        if _VERBOSE:
            print("    get_gtests", gtestfile)
        with open(gtestfile) as file:
            contents = file.read()
            for match in re.findall(
                    "TEST\(\s*(\S+?),\s*(\S+?)\s*\)", contents):
                assert len(match) is 2
                gtests.add(match[0])
    return gtests


def src_dir_new_file():
    fname = '.__build_utils_src_dir_contents'
    ftemp = '.__build_utils_src_dir_contents_tmp'
    os.system('find ./src -name \\*.cpp > ' + ftemp)
    are_same = filecmp.cmp(fname, ftemp) if exists(fname) else False
    if not are_same:
        os.rename(ftemp, fname)
    else:
        os.remove(ftemp)
    return not are_same


def filter_testfiles(testfiles):
    testfiles = [x for x in testfiles if '/tools/' not in x]
    return testfiles


def rebuild_fast(target='_rif', cfg='Release', force_redo_cmake=False):
    try:
        if force_redo_cmake:
            raise OSError
        cmake_dir = get_cmake_dir('temp', cfg=cfg)
        makeexe = 'ninja'
        if not which('ninja'):
            makeexe = 'make'
        ncpu = multiprocessing.cpu_count()
        return os.system('cd ' + cmake_dir + '; ' +
                         makeexe + ' -j%i ' % ncpu + target)
    except (OSError, TypeError):
        return rebuild_rif(cfg=cfg)


def build_and_test():
    assert exists('CMakeLists.txt') and exists('setup.py')
    # print("== build_and_test ==")
    argfiles = sys.argv[1:]
    cfg = 'Release'
    srcdir = abspath('src')
    testfiles = set()
    for argfile in argfiles:
        if not argfile.endswith('.py') and not argfile.startswith('-'):
            # testfiles.add(argfile)
            continue
        candidate_mod = argfile.replace('_test.py', '.py')
        candidate_tst = candidate_mod.replace('.py', '_test.py')
        if exists(candidate_mod) and exists(candidate_tst):
            testfiles.add(candidate_mod)
            testfiles.add(candidate_tst)
    testfiles = list(testfiles)
    # print(testfiles)
    # sys.exit(-1)
    testfiles.extend(x.replace('.pybind.cpp', '_test.py') for x in argfiles
                     if x.endswith('.pybind.cpp') and
                     exists(x.replace('.pybind.cpp', '_test.py')))
    testfiles.extend(x.replace('.pybind.cpp', '.py') for x in argfiles
                     if x.endswith('.pybind.cpp') and
                     exists(x.replace('.pybind.cpp', '.py')))
    for x in argfiles:
        if x.endswith('.hpp'):
            if not exists(x.replace('.hpp', '.gtest.cpp')):
                candidate = x.replace('.hpp', '_test.py')
                if exists(candidate):
                    testfiles.append(candidate)
    pybindfiles = [x for x in argfiles if x.endswith('.pybind.cpp')]
    gtests = get_gtests(argfiles)
    apps = [x for x in argfiles if '/rif/apps/' in x]
    tests = [x for x in argfiles if '/rif/apps/' not in x]
    # print('== calling rebuild_fast ==')
    no_xdist = len(testfiles) or len(gtests)
    no_xdist |= os.system('egrep "#.*cpp_files" pytest.ini') == 0
    # no_xdist |= sys.version_info.major is 3 and sys.version_info.minor is 6
    force_redo_cmake = len(pybindfiles) or not no_xdist or src_dir_new_file()
    rebuild_fast(target='_rif gtest_all',
                 cfg=cfg, force_redo_cmake=force_redo_cmake)
    libdir = abspath(get_cmake_dir('lib', cfg))
    builddir = dirname(libdir)
    testfiles = filter_testfiles(testfiles)
    testfiles_in_lib = [f.replace(srcdir, libdir) for f in testfiles]
    for a, b in zip(testfiles, testfiles_in_lib):
        shutil.copyfile(a, b)
        if exists(a.replace('_test.py', '.py')):
            shutil.copyfile(a.replace('_test.py', '.py'),
                            b.replace('_test.py', '.py'))
    shutil.copy(join(srcdir, 'rif', 'conftest.py'),
                join(libdir, 'rif', 'conftest.py'))
    testfiles = testfiles_in_lib

    if '--noxdist' in sys.argv:
        no_xdist = True

    if '--inplace' in sys.argv:
        if _VERBOSE:
            print('== adding to python path:', libdir, '==')
        assert exists(libdir)
        # need to use sys.path for this process
        sys.path.insert(0, libdir)
        # need to use PYTHONPATH env for xdist subprocesses
        add_to_pypath(libdir)
        assert libdir in sys.path
        assert libdir in os.environ['PYTHONPATH'].split(':')

    if apps:
        for app in apps:
            sys.stdout.write('===== running:' + app + '=====' + os.linesep)
            sys.stdout.flush()
            os.system(sys.executable + ' ' + app)
            sys.stdout.write('===== done: ' + app + ' ======' + os.linesep)
            sys.stdout.flush()
    if tests or not apps:
        ncpu = get_ncpu()
        args = list(testfiles)
        if not (testfiles or gtests):
            args = ['buildR', '--pyargs', 'rif', ]
        args.extend(['--doctest-modules'])
        args.extend(['--ignore', 'build'])
        args.extend(['--ignore', 'src'])
        if not no_xdist:
            # args.extend('--cov=./src -n{}'.format(ncpu).split())
            args.extend('-n{}'.format(min(ncpu, 16)).split())
        if gtests:
            args.extend(['-k', ' or '.join(gtests)])
        if testfiles and not gtests:
            args.extend(['--ignore', builddir])
        for decoy in get_ignored_dirs(cfg):
            args += ['--ignore', decoy]
        if _VERBOSE:
            print('==============================================================')
            sys.stdout.write('pytest' + os.linesep + '   ')
            for arg in args:
                if arg.startswith('-'):
                    sys.stdout.write(os.linesep + '   ')
                sys.stdout.write(' ' + arg)
            print()
            print('    sys.path:')
            for p in sys.path:
                print('       ', p)
            print('    cwd:', os.getcwd())
        sys.argv[1:] = args
        assert not pytest.main()


def make_gtest_auto_cpp(files, cmake_dir):
    # print('make_gtest_auto_cpp!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
    includes = os.linesep.join("#include <%s>" % f for f in files)
    code = """#include "test/gtest_util.hpp"
%s
int main(int argc, char **argv) {
  std::vector<std::string> args;
  for (int i = 0; i < argc; ++i) args.push_back(std::string(argv[i]));
  // std::cout << "int main(int argc, char **argv) FROM " << __FILE__ << std::endl;
  init_gtest_tests(args);
  return run_gtest_tests();
}
""" % includes
    if includes:
        if _VERBOSE:
            print('build_and_run_gtest_auto.py: making gtest_auto.gen.cpp')
        with open(cmake_dir + '/gtest_auto.gen.cpp', 'w') as out:
            out.write(code)


def build_and_run_gtest_auto():
    try:
        cmake_dir = get_cmake_dir('temp', cfg='Release')
        assert cmake_dir
    except AssertionError:
        rebuild_fast('gtest_wip')
        cmake_dir = get_cmake_dir('temp', cfg='Release')
    files = [x for x in sys.argv[1:] if x.endswith('.gtest.cpp')]
    for x in sys.argv[1:]:
        if x.endswith('.hpp'):
            x = x.replace('.hpp', '.gtest.cpp')
            candidates = [x]
            for i in range(x.count('_')):
                candidates.append(
                    '_'.join(x.split('_')[:i + 1]) + '.gtest.cpp')
            for f in candidates:
                if exists(f):
                    files.append(f)
    files = set(files)
    # print('    files:', files)
    # nargs = len([x for x in sys.argv if not x.startswith('-')])
    if not len(files):  # or len(files) + 1 != nargs:
        raise NotImplementedError
    make_gtest_auto_cpp(files, cmake_dir)
    assert not os.system('cd ' + cmake_dir + ' && ninja gtest_auto')
    assert not os.system(cmake_dir + '/gtest_auto')
