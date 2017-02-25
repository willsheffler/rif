from builtins import bytes

import glob
import os
import sys


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


def get_build_dir(d, cfg='Release'):
    """get directory setup.py builds stuff in, d is lib or temp"""
    version = '{}.{}'.format(sys.version_info.major, sys.version_info.minor)
    libdir = (glob.glob(get_proj_root() + '/build_setup_py_' +
                        cfg + '/' + d + '*' + version))
    assert len(libdir) < 2
    if libdir:
        return libdir[0]
    else:
        return None


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
    cmake_dir = get_build_dir('temp', cfg=cfg)
    if not cmake_dir or redo_cmake:
        if rebuild_setup_py_rif(cfg=cfg):
            return -1
        cmake_dir = get_build_dir('temp', cfg=cfg)

    return os.system('cd ' + cmake_dir + '; ' + makeexe + ' -j8 ' + target)


def make_docs(kind='html'):
    proj_root = get_proj_root()
    rebuild_setup_py_rif()
    add_to_pypath(get_build_dir('lib'))
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
    # remove_installed_rif()
    proj_root = get_proj_root()
    print('calling rebuild_fast')
    if rebuild_fast(target='rif_cpp gtest_all',
                    cfg='Release', redo_cmake=redo_cmake):
        sys.exit(-1)
    pypath = [os.path.abspath(get_build_dir('lib')),
              proj_root + '/external',
              ]
    # need to use sys.path for this process
    sys.path = pypath + sys.path
    # need to use PYTHONPATH env for xdist subprocessess
    add_to_pypath(pypath)
    # TODO both here and in docs, this gets messed
    #      up when rif is actually installed
    import pytest
    if sys.version_info.major is 2:
        proj_root = bytes(proj_root, 'ascii')
    args = [x for x in sys.argv[1:] if x.endswith('.py') and
            os.path.basename(x).startswith('test')]
    if not args:
        args = ['.', '-n4']
    pytest.main(args)
