from rif.vis import showme
import pytest
from rif import rcl

try:
    import pymol
    HAVE_PYMOL = True
except ImportError:
    HAVE_PYMOL = False


def test_showme_error():
    with pytest.raises(NotImplementedError):
        showme('foo', how='text')


@pytest.mark.skipif('not rcl.HAVE_PYROSETTA or not HAVE_PYMOL')
def test_showme_rosetta_pymol(curved_helix_pose):
    pose = curved_helix_pose
    things = [pose] * 10
    # showme(things, headless=1)
