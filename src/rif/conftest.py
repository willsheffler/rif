import pytest
import os
from os.path import join, dirname, abspath, exists
from rif import rcl


@pytest.fixture(scope='session')
def rootdir():
    d = abspath(dirname(__file__))
    assert exists(d)
    return d


@pytest.fixture(scope='session')
def datadir(rootdir):
    d = join(rootdir, 'data')
    assert exists(d)
    return d


@pytest.fixture(scope='session')
def pdbdir(datadir):
    d = join(datadir, 'pdb')
    assert exists(d)
    return d


@pytest.fixture(scope='session')
def pdbfname(datadir):
    f = join(datadir, 'pdb', '3uc7A.pdb')
    assert exists(f)
    return f


@pytest.fixture(scope='session')
def cleanpdbfname(pdbdir):
    f = join(pdbdir, '3uc7A_clean.pdb')
    assert exists(f)
    return f


@pytest.fixture(scope='session')
def pdbsmallfname(pdbdir):
    f = join(pdbdir, 'small.pdb')
    assert exists(f)
    return f


@pytest.fixture(scope='session')
def pose(pdbsmallfname):
    if not rcl.HAVE_PYROSETTA: return None
    rcl.init_check(strict=False)
    pose = rcl.pose_from_file(pdbsmallfname)
    return pose


@pytest.fixture(scope='session')
def bigpose(cleanpdbfname):
    if not rcl.HAVE_PYROSETTA: return None
    rcl.init_check(strict=False)
    pose = rcl.pose_from_file(cleanpdbfname)
    return pose


@pytest.fixture(scope='session')
def curved_helix_pose(pdbdir):
    if not rcl.HAVE_PYROSETTA: return None
    rcl.init_check(strict=False)
    pose = rcl.pose_from_file(join(pdbdir, 'curved_helix.pdb'))
    return pose


@pytest.fixture(scope='session')
def strand_pose(pdbdir):
    if not rcl.HAVE_PYROSETTA: return None
    rcl.init_check(strict=False)
    pose = rcl.pose_from_file(join(pdbdir, 'strand.pdb'))
    return pose


@pytest.fixture(scope='session')
def loop_pose(pdbdir):
    if not rcl.HAVE_PYROSETTA: return None
    rcl.init_check(strict=False)
    pose = rcl.pose_from_file(join(pdbdir, 'loop.pdb'))
    return pose


def get_smalltrimer(pdbdir, suffix=''):
    if not rcl.HAVE_PYROSETTA: return None
    rcl.init_check(strict=False)
    return rcl.pose_from_file(join(pdbdir, '1coi%s.pdb' % suffix))


@pytest.fixture(scope='session')
def small_trimer_A(pdbdir):
    return get_smalltrimer(pdbdir, suffix='_A')


@pytest.fixture(scope='session')
def small_trimer_B(pdbdir):
    return get_smalltrimer(pdbdir, suffix='_B')


@pytest.fixture(scope='session')
def small_trimer_C(pdbdir):
    return get_smalltrimer(pdbdir, suffix='_C')
