import pytest
from rif import rcl
from rif.data import poselib as plib


@pytest.mark.skipif('not rcl.HAVE_PYROSETTA')
def test_poselib():
    assert len(plib.curved_helix) is 13
    assert len(plib.strand) is 9
    assert len(plib.loop) is 8
    assert len(plib.small) is 7
    assert len(plib.get('1coi_A')) is 29
    assert plib.strand is plib.strand  # check cache...
