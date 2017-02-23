import pytest
from rif.visual.vispymol import launch_pymol


@pytest.mark.skip
def test_launch_pymol():
    def foo():
        pass
    launch_pymol(foo, args='-qc')
