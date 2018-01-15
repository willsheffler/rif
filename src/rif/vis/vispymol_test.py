from rif.vis import showme, pymol_xform
import pytest
from rif import rcl
import numpy as np
try:
    import pymol
    HAVE_PYMOL = True
except ImportError:
    HAVE_PYMOL = False


def test_showme_error():
    with pytest.raises(NotImplementedError):
        showme('foo', how='text')


# @pytest.mark.skipif('not rcl.HAVE_PYROSETTA or not HAVE_PYMOL')
# def test_showme_rosetta_pymol(curved_helix_pose):
#     pose = curved_helix_pose
#     things = [pose] * 10
#     showme(things, headless=1)


# @pytest.mark.skipif('not rcl.HAVE_PYROSETTA or not HAVE_PYMOL')
# def test_showme_rosetta_pymol(curved_helix_pose):
#     pose = curved_helix_pose
#     name = showme(pose, headless=1)['last_obj']
#     x = np.array([[1, 0, 0, 100],
#                   [0, 1, 0, 200],
#                   [0, 0, 1, 300],
#                   [0, 0, 0, 1]])
#     import pymol
#     xyz0 = pymol.com('all')
#     pymol_xform(name, x)
#     xyz = pymol.com('all')
#     print(xyz0)
#     print(xyz)
#     assert abs((xyz - xyz0).x - 100) < 0.001
#     assert abs((xyz - xyz0).y - 200) < 0.001
#     assert abs((xyz - xyz0).z - 300) < 0.001
