from __future__ import absolute_import, division, print_function
from builtins import bytes

import glob
import os
import sys
import multiprocessing

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
    if len(libdir) > 1:
        print('ERROR get_cmake_dir', path)
        print(libdir)
        assert len(libdir) < 2
    if libdir:
        return libdir[0]
    else:
        print('cant find libdir', path)
        raise IOError


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
    if os.system('cd ' + proj_root + '; ' + sys.executable +
                 ' setup.py build --build-base=build_setup_py_' + cfg):
        return -1


def rebuild_fast(target='rif_cpp', cfg='Release', redo_cmake=False):
    makeexe = 'ninja'
    if not which('ninja'):
        makeexe = 'make'
    # proj_root = get_proj_root()
    try:
        cmake_dir = get_cmake_dir('temp', cfg=cfg)
    except:
        cmake_dir = None
    if not cmake_dir or redo_cmake:
        if rebuild_setup_py_rif(cfg=cfg):
            return -1
        cmake_dir = get_cmake_dir('temp', cfg=cfg)

    return os.system('cd ' + cmake_dir + '; ' + makeexe + ' -j8 ' + target)


def make_docs(kind='html', cfg='Release'):
    proj_root = get_proj_root()
    rebuild_setup_py_rif()
    add_to_pypath(get_cmake_dir('lib', cfg=cfg))
    os.system('cd ' + proj_root + '/docs; make ' + kind)


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


def rif_is_installed():
    raise NotImplemented


def remove_installed_rif():
    if rif_is_installed():
        os.system('echo y | ' + sys.executable + ' -m pip uninstall rif')
        assert is_devel_install() or not rif_is_installed()
        print('uninstalled rif from ' + sys.executable)


def build_and_run_pytest(redo_cmake=False):
    cfg = 'Release'
    # remove_installed_rif()
    proj_root = get_proj_root()
    build_dir = get_build_dir(cfg)
    print('calling rebuild_fast')
    if not os.path.exists(build_dir):
        rebuild_setup_py_rif(cfg)
    assert os.path.exists(build_dir)
    if rebuild_fast(target='rif_cpp gtest_all',
                    cfg=cfg, redo_cmake=redo_cmake):
        sys.exit(-1)
    # TODO both here and in docs, this gets messed
    #      up when rif is actually installed
    libdir = os.path.abspath(get_cmake_dir('lib', cfg))
    print('adding to python path:', libdir)
    assert os.path.exists(libdir)
    # need to use sys.path for this process
    sys.path.append(libdir)
    # need to use PYTHONPATH env for xdist subprocessess
    add_to_pypath(libdir)
    assert libdir in sys.path
    assert libdir in os.environ['PYTHONPATH'].split(':')
    if sys.version_info.major is 2:
        proj_root = bytes(proj_root, 'ascii')
    args = [x for x in sys.argv[1:] if x.endswith('.py') and
            os.path.basename(x).startswith('test')]
    ncpu = multiprocessing.cpu_count()
    if ncpu > 2:
        ncpu = int(ncpu / 2)
    if 'CI' in os.environ:
        ncpu = 4
    if not args:
        args = '. --ignore build --cov=./src -n{}'.format(ncpu).split()
    else:  # running one file, don't scan
        args += '--ignore build_setup_py_Release'.split()
    for decoy in get_ignored_dirs(cfg):
        args += ['--ignore', decoy]
    print('pytest.main(', ' '.join(args), ')')
    return pytest.main(args)
