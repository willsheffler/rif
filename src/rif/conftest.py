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


def get_pose(pdbdir, fname):
    if not rcl.HAVE_PYROSETTA: return None
    rcl.init_check(strict=False)
    pose = rcl.pose_from_file(join(pdbdir, fname))
    return pose


@pytest.fixture(scope='session')
def curved_helix_pose(pdbdir):
    return get_pose(pdbdir, 'curved_helix.pdb')


@pytest.fixture(scope='session')
def strand_pose(pdbdir):
    return get_pose(pdbdir, 'strand.pdb')


@pytest.fixture(scope='session')
def loop_pose(pdbdir):
    return get_pose(pdbdir, 'loop.pdb')


@pytest.fixture(scope='session')
def trimer_pose(pdbdir):
    return get_pose(pdbdir, '1coi.pdb')


@pytest.fixture(scope='session')
def trimerA_pose(pdbdir):
    return get_pose(pdbdir, '1coi_A.pdb')


@pytest.fixture(scope='session')
def trimerB_pose(pdbdir):
    return get_pose(pdbdir, '1coi_B.pdb')


@pytest.fixture(scope='session')
def trimerC_pose(pdbdir):
    return get_pose(pdbdir, '1coi_C.pdb')


@pytest.fixture(scope='session')
def c2pose(pdbdir):
    return get_pose(pdbdir, 'C2_4agh_1_full.pdb')


@pytest.fixture(scope='session')
def c3pose(pdbdir):
    return get_pose(pdbdir, 'C3_1wp8_1_full.pdb')


@pytest.fixture(scope='session')
def c4pose(pdbdir):
    return get_pose(pdbdir, 'C4_1gcl_1_full.pdb')


@pytest.fixture(scope='session')
def c5pose(pdbdir):
    return get_pose(pdbdir, 'C5_3mxg_1_full.pdb')


@pytest.fixture(scope='session')
def c6pose(pdbdir):
    return get_pose(pdbdir, 'C6_2xf5_1_full.pdb')
