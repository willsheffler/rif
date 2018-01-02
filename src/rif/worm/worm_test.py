from rif import rcl
from rif.homog import hrot, htrans, axis_angle_of
from rif.vis import showme
from numpy.testing import assert_allclose
import pytest
import numpy as np
if rcl.HAVE_PYROSETTA:
    from rif.worm import *
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor


@pytest.mark.skipif('not rcl.HAVE_PYROSETTA')
def test_SpliceSite(pose):
    assert len(pose) == 7
    ss = SpliceSite(1, 'N')
    assert 1 == ss.resid(1, pose)
    assert pose.size() == ss.resid(-1, pose)
    assert ss.resids(pose) == [1]
    assert SpliceSite('1:7', 'N').resids(pose) == [1, 2, 3, 4, 5, 6, 7]
    assert SpliceSite(':7', 'N').resids(pose) == [1, 2, 3, 4, 5, 6, 7]
    assert SpliceSite('-3:-1', 'N').resids(pose) == [5, 6, 7]
    assert SpliceSite('-3:', 'N').resids(pose) == [5, 6, 7]
    assert SpliceSite(':2', 'N').resids(pose) == [1, 2]
    assert SpliceSite(':-5', 'N').resids(pose) == [1, 2, 3]
    assert SpliceSite('::2', 'N').resids(pose) == [1, 3, 5, 7]
    with pytest.raises(ValueError): SpliceSite('-1:-3', 'N').resids(pose)
    with pytest.raises(ValueError): SpliceSite('-1:3', 'N').resids(pose)


def test_geom_check():
    SX = SegmentSym
    I = np.identity(4, 'f4')
    rotx1rad = hrot([1, 0, 0], 1)
    transx10 = htrans([10, 0, 0])
    randaxes = np.random.randn(1, 3)

    assert 0 == SX('c1').score([I, I])
    assert 0.001 > abs(100 - SX('c1').score([I, rotx1rad]))
    assert 1e-5 > abs(SX('c2').score([I, hrot([1, 0, 0], np.pi)]))

    score = SegmentSym('c2').score([I, hrot(randaxes, np.pi)])
    assert_allclose(0, score, atol=1e-5, rtol=1)

    score = SegmentSym('c3').score([I, hrot(randaxes, np.pi * 2 / 3)])
    assert_allclose(0, score, atol=1e-5, rtol=1)

    score = SegmentSym('c4').score([I, hrot(randaxes, np.pi / 2)])
    assert_allclose(0, score, atol=1e-5, rtol=1)


@pytest.mark.skipif('not rcl.HAVE_PYROSETTA')
def test_segment_geom(curved_helix_pose):
    "currently only a basic sanity checkb... only checks translation distances"
    body = curved_helix_pose
    stubs = rcl.bbstubs(body)['raw']
    assert stubs.shape == (body.size(), 4, 4)

    nsplice = SpliceSite(polarity='N', sele=[1, 2, ])
    csplice = SpliceSite(polarity='C', sele=[9, 10, 11, 12, 13])
    Npairs0 = len(nsplice.selections) * len(csplice.selections)

    # N to N and C to C invalid, can't splice to same
    splicable = Spliceable(body, sites=[nsplice, csplice])
    with pytest.raises(ValueError):
        seg = Segment([splicable], entry='N', exit='N')
    with pytest.raises(ValueError):
        seg = Segment([splicable] * 3, entry='C', exit='C')

    # add some extra splice sites
    Nexsite = 2
    splicable = Spliceable(body, sites=[nsplice, csplice] * Nexsite)

    # test beginning segment.. only has exit
    seg = Segment([splicable], exit='C')
    assert seg.x2exit.shape == (Nexsite * len(csplice.selections), 4, 4)
    assert seg.x2orgn.shape == (Nexsite * len(csplice.selections), 4, 4)
    assert np.all(seg.x2exit[..., 3, :3] == 0)
    assert np.all(seg.x2exit[..., 3, 3] == 1)
    for e2x, e2o, ir, jr in zip(seg.x2exit, seg.x2orgn,
                                seg.entryresid, seg.exitresid):
        assert ir == -1
        assert_allclose(e2o, np.eye(4))
        assert_allclose(e2x, stubs[jr - 1])

    # test middle segment with entry and exit
    seg = Segment([splicable], entry='N', exit='C')
    assert seg.x2exit.shape == (Nexsite**2 * Npairs0, 4, 4)
    assert seg.x2orgn.shape == (Nexsite**2 * Npairs0, 4, 4)
    assert np.all(seg.x2exit[..., 3, :3] == 0)
    assert np.all(seg.x2exit[..., 3, 3] == 1)
    for e2x, e2o, ir, jr in zip(seg.x2exit, seg.x2orgn,
                                seg.entryresid, seg.exitresid):
        assert_allclose(stubs[ir - 1] @ e2o, np.eye(4), atol=1e-5)
        assert_allclose(stubs[ir - 1] @ e2x, stubs[jr - 1], atol=1e-5)

    # test ending segment.. only has entry
    seg = Segment([splicable], entry='N')
    assert seg.x2exit.shape == (Nexsite * len(nsplice.selections), 4, 4)
    assert seg.x2orgn.shape == (Nexsite * len(nsplice.selections), 4, 4)
    assert np.all(seg.x2exit[..., 3, :3] == 0)
    assert np.all(seg.x2exit[..., 3, 3] == 1)
    for e2x, e2o, ir, jr in zip(seg.x2exit, seg.x2orgn,
                                seg.entryresid, seg.exitresid):
        assert jr == -1
        assert_allclose(e2o, e2x)
        assert_allclose(e2o @ stubs[ir - 1], np.eye(4), atol=1e-5)

    # test now with multiple splicables input to segment
    Nexbody = 3
    seg = Segment([splicable] * Nexbody, entry='N', exit='C')
    Npairs_expected = Nexbody * Nexsite**2 * Npairs0
    assert seg.x2exit.shape == (Npairs_expected, 4, 4)
    assert seg.x2orgn.shape == (Npairs_expected, 4, 4)
    assert len(seg.entryresid) == Npairs_expected
    assert len(seg.exitresid) == Npairs_expected
    assert len(seg.bodyid) == Npairs_expected
    for i in range(Nexbody):
        assert i == seg.bodyid[0 + i * Npairs0 * Nexsite**2]
    assert np.all(seg.x2exit[..., 3, :3] == 0)
    assert np.all(seg.x2exit[..., 3, 3] == 1)
    for e2x, e2o, ir, jr in zip(seg.x2exit, seg.x2orgn,
                                seg.entryresid, seg.exitresid):
        assert_allclose(stubs[ir - 1] @ e2o, np.eye(4), atol=1e-5)
        assert_allclose(stubs[ir - 1] @ e2x, stubs[jr - 1], atol=1e-5)

        # this is incorrect... translation includes
        # cart motion due to non-zero rotation center
        # ijdis = np.sqrt(np.sum((e2x[..., :3, 3]**2)))
        # cadis = body.residue(ir).xyz('CA').distance(body.residue(jr).xyz('CA'))
        # assert abs(ijdis - cadis) < 0.001
        # odis = np.sqrt(np.sum((e2o[..., :3, 3]**2)))
        # caodis = body.residue(ir).xyz('CA').length()
        # assert abs(odis - caodis) < 0.001


@pytest.mark.skipif('not rcl.HAVE_PYROSETTA')
def test_grow_cycle(curved_helix_pose):
    helix = Spliceable(curved_helix_pose, sites=[(1, 'N'), ('-4:', 'C')])
    segments = ([Segment([helix], exit='C'), ] +
                [Segment([helix], entry='N', exit='C')] * 3 +
                [Segment([helix], entry='N')])
    worms = grow(segments, SegmentSym('C2', lever=20))
    assert 0.1411 < np.min(worms.scores) < 0.1412


@pytest.mark.skipif('not rcl.HAVE_PYROSETTA')
def test_grow_cycle_thread_pool(curved_helix_pose):
    helix = Spliceable(curved_helix_pose, sites=[(1, 'N'), ('-4:', 'C')])
    segments = ([Segment([helix], exit='C'), ] +
                [Segment([helix], entry='N', exit='C')] * 3 +
                [Segment([helix], entry='N')])
    worms = grow(segments, SegmentSym('C2', lever=20),
                 executor=ThreadPoolExecutor, max_workers=2)
    assert 0.1411 < np.min(worms.scores) < 0.1412
    assert np.sum(worms.scores < 0.1412) == 4


@pytest.mark.skipif('not rcl.HAVE_PYROSETTA')
def test_grow_cycle_process_pool(curved_helix_pose):
    helix = Spliceable(curved_helix_pose, sites=[(1, 'N'), ('-4:', 'C')])
    segments = ([Segment([helix], exit='C'), ] +
                [Segment([helix], entry='N', exit='C')] * 3 +
                [Segment([helix], entry='N')])
    worms = grow(segments, SegmentSym('C2', lever=20),
                 executor=ProcessPoolExecutor, max_workers=2)
    assert 0.1411 < np.min(worms.scores) < 0.1412
    assert np.sum(worms.scores < 0.1412) == 4


@pytest.mark.skipif('not rcl.HAVE_PYROSETTA')
def test_grow_errors(curved_helix_pose):
    nsplice = SpliceSite(sele=[1, 2, 3, 4, 5, 6], polarity='N')
    csplice = SpliceSite(sele=[13, ], polarity='C')
    splicable1 = Spliceable(body=curved_helix_pose, sites=[nsplice, csplice])
    splicable2 = Spliceable(body=curved_helix_pose, sites=[nsplice, csplice])
    splicables = [splicable1]
    segments = ([Segment(splicables, exit='C'), ] +
                [Segment(splicables, entry='N', exit='C'), ] * 3 +
                [Segment(splicables, entry='N'), ])
    checkc3 = SegmentSym('C2', from_seg=0, to_seg=-1)

    # make sure incorrect begin/end throws error
    with pytest.raises(ValueError):
        grow(segments[: 2], criteria=checkc3)
    with pytest.raises(ValueError):
        grow(segments[1:], criteria=checkc3)
    segments_polarity_mismatch = [
        Segment(splicables, exit='C'),
        Segment(splicables, entry='C'),
    ]
    with pytest.raises(ValueError):
        grow(segments_polarity_mismatch, criteria=checkc3)


@pytest.mark.skipif('not rcl.HAVE_PYROSETTA')
def test_memlimit(curved_helix_pose):
    helix = Spliceable(curved_helix_pose, sites=[((1, 2), 'N'), ('-2:', 'C')])
    segments = ([Segment([helix], exit='C'), ] +
                [Segment([helix], entry='N', exit='C')] * 3 +
                [Segment([helix], entry='N')])
    for i in range(2, 7):
        w1 = grow(segments, SegmentSym('c2'), memlim=10**i, thresh=30)
        assert i == 2 or np.allclose(w0.scores, w1.scores)
        w0 = w1


@pytest.mark.skipif('not rcl.HAVE_PYROSETTA')
def test_pose_alignment(curved_helix_pose):
    helix = Spliceable(curved_helix_pose, sites=[(1, 'N'), ('-4:', 'C')])
    segments = ([Segment([helix], exit='C'), ] +
                [Segment([helix], entry='N', exit='C')] * 3 +
                [Segment([helix], entry='N')])
    w = grow(segments, SegmentSym('c2'), thresh=1)
    assert len(w)
    pose = w.pose(0, onechain=True, align=1, withend=1)
    xyz0 = np.array([pose.residue(1).xyz(2)[i] for i in (0, 1, 2)] + [1])
    xyz1 = np.array([pose.residue(43).xyz(2)[i] for i in (0, 1, 2)] + [1])
    xyz1 = hrot([0, 0, 1], 180) @ xyz1
    assert np.sum((xyz1 - xyz0)**2) < 0.1
