from rif.worm import *
from rif.vis import showme
from rif import rcl
import pytest

def test_target_geometry():
    assert 1


@pytest.mark.skipif('not rcl.HAVE_PYROSETTA')
def test_worm(curved_helix_pose):
    # showme(curved_helix_pose)
    assert 1
