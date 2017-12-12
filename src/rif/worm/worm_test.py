from rif.worm import *
from rif.vis import showme
from rif import rcl
import pytest


def test_target_geometry():
    assert 1


@pytest.mark.skipif('not rcl.HAVE_PYROSETTA')
def test_segment_geom(curved_helix_pose):
    pose = curved_helix_pose
    nsplice = SplicePositions([1, 2, 3, 4], 'N')
    csplice = SplicePositions([9, 10, 11, 12, 13], 'C')
    splice_groups = (nsplice, csplice)
    with pytest.raises(AssertionError):
        # no valid splices here, only one N can't connect to itself
        s = Segment(pose, splice_groups, 'N', 'N')

    # test beginning segment.. only has exit
    s = Segment(pose, splice_groups, exit_dir='C')
    assert s.entry2exit.shape == (len(csplice), 4, 4)
    assert (s.entry2exit[..., 3, :3] == 0).all()
    assert (s.entry2exit[..., 3, 3] == 1).all()
    for x, ir, jr in zip(s.entry2exit, s.entry_resi, s.exit_resi):
        assert np.isnan(ir)
        jdis = np.sqrt(np.sum((x[..., :3, 3]**2)))
        cadis = pose.residue(jr).xyz('CA').length()
        assert abs(jdis - cadis) < 0.001

    # test middle segment with entry and exit
    s = Segment(pose, splice_groups, entry_dir='N', exit_dir='C')
    assert s.entry2exit.shape == (len(nsplice) * len(csplice), 4, 4)
    assert (s.entry2exit[..., 3, :3] == 0).all()
    assert (s.entry2exit[..., 3, 3] == 1).all()
    for x, ir, jr in zip(s.entry2exit, s.entry_resi, s.exit_resi):
        ijdis = np.sqrt(np.sum((x[..., :3, 3]**2)))
        cadis = pose.residue(ir).xyz('CA').distance(pose.residue(jr).xyz('CA'))
        assert abs(ijdis - cadis) < 0.001

    # test ending segment.. only has entry
    s = Segment(pose, splice_groups, entry_dir='N')
    assert s.entry2exit.shape == (len(nsplice), 4, 4)
    assert (s.entry2exit[..., 3, :3] == 0).all()
    assert (s.entry2exit[..., 3, 3] == 1).all()
    for x, ir, jr in zip(s.entry2exit, s.entry_resi, s.exit_resi):
        assert np.isnan(jr)
        idis = np.sqrt(np.sum((x[..., :3, 3]**2)))
        cadis = pose.residue(ir).xyz('CA').length()
        assert abs(idis - cadis) < 0.001
