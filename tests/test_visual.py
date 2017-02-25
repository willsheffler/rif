import pytest
from rif.visual.vispymol import launch_pymol, HAVE_PYMOL


@pytest.mark.skip
@pytest.mark.skipif('not HAVE_PYMOL')
def test_launch_pymol():
    def foo():
        pass
    launch_pymol(foo, args='-qc')
