from rif import rcl
from rif.homog import hrot, htrans, axis_angle_of, axis_ang_cen_of
from rif.vis.vispymol import showme, showline, showsphere
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
    spliceable = Spliceable(pose, [])
    assert 1 == ss.resid(1, spliceable.body)
    assert pose.size() == ss.resid(-1, spliceable.body)
    assert ss.resids(spliceable) == [1]
    assert SpliceSite('1:7', 'N').resids(spliceable) == [1, 2, 3, 4, 5, 6, 7]
    assert SpliceSite(':7', 'N').resids(spliceable) == [1, 2, 3, 4, 5, 6, 7]
    assert SpliceSite('-3:-1', 'N').resids(spliceable) == [5, 6, 7]
    assert SpliceSite('-3:', 'N').resids(spliceable) == [5, 6, 7]
    assert SpliceSite(':2', 'N').resids(spliceable) == [1, 2]
    assert SpliceSite(':-5', 'N').resids(spliceable) == [1, 2, 3]
    assert SpliceSite('::2', 'N').resids(spliceable) == [1, 3, 5, 7]
    with pytest.raises(ValueError): SpliceSite('-1:-3', 'N').resids(spliceable)
    with pytest.raises(ValueError): SpliceSite('-1:3', 'N').resids(spliceable)


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
    spliceable = Spliceable(body, sites=[nsplice, csplice])
    with pytest.raises(ValueError):
        seg = Segment([spliceable], entry='N', exit='N')
    with pytest.raises(ValueError):
        seg = Segment([spliceable] * 3, entry='C', exit='C')

    # add some extra splice sites
    Nexsite = 2
    spliceable = Spliceable(body, sites=[nsplice, csplice] * Nexsite)

    # test beginning segment.. only has exit
    seg = Segment([spliceable], exit='C')
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
    seg = Segment([spliceable], 'N', 'C')
    assert seg.x2exit.shape == (Nexsite**2 * Npairs0, 4, 4)
    assert seg.x2orgn.shape == (Nexsite**2 * Npairs0, 4, 4)
    assert np.all(seg.x2exit[..., 3, :3] == 0)
    assert np.all(seg.x2exit[..., 3, 3] == 1)
    for e2x, e2o, ir, jr in zip(seg.x2exit, seg.x2orgn,
                                seg.entryresid, seg.exitresid):
        assert_allclose(stubs[ir - 1] @ e2o, np.eye(4), atol=1e-5)
        assert_allclose(stubs[ir - 1] @ e2x, stubs[jr - 1], atol=1e-5)

    # test ending segment.. only has entry
    seg = Segment([spliceable], entry='N')
    assert seg.x2exit.shape == (Nexsite * len(nsplice.selections), 4, 4)
    assert seg.x2orgn.shape == (Nexsite * len(nsplice.selections), 4, 4)
    assert np.all(seg.x2exit[..., 3, :3] == 0)
    assert np.all(seg.x2exit[..., 3, 3] == 1)
    for e2x, e2o, ir, jr in zip(seg.x2exit, seg.x2orgn,
                                seg.entryresid, seg.exitresid):
        assert jr == -1
        assert_allclose(e2o, e2x)
        assert_allclose(e2o @ stubs[ir - 1], np.eye(4), atol=1e-5)

    # test now with multiple spliceables input to segment
    Nexbody = 3
    seg = Segment([spliceable] * Nexbody, 'N', 'C')
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


@pytest.mark.skipif('not rcl.HAVE_PYROSETTA')
def test_grow_cycle(curved_helix_pose):
    helix = Spliceable(curved_helix_pose, sites=[(1, 'N'), ('-4:', 'C')])
    segments = ([Segment([helix], exit='C'), ] +
                [Segment([helix], 'N', 'C')] * 3 +
                [Segment([helix], entry='N')])
    worms = grow(segments, SegmentSym('C2', lever=20))
    assert 0.1411 < np.min(worms.scores) < 0.1412


@pytest.mark.skipif('not rcl.HAVE_PYROSETTA')
def test_grow_cycle_thread_pool(curved_helix_pose):
    helix = Spliceable(curved_helix_pose, sites=[(1, 'N'), ('-4:', 'C')])
    segments = ([Segment([helix], exit='C'), ] +
                [Segment([helix], 'N', 'C')] * 3 +
                [Segment([helix], entry='N')])
    worms = grow(segments, SegmentSym('C2', lever=20),
                 executor=ThreadPoolExecutor, max_workers=2)
    assert 0.1411 < np.min(worms.scores) < 0.1412
    assert np.sum(worms.scores < 0.1412) == 4


@pytest.mark.skipif('not rcl.HAVE_PYROSETTA')
def test_grow_cycle_process_pool(curved_helix_pose):
    helix = Spliceable(curved_helix_pose, sites=[(1, 'N'), ('-4:', 'C')])
    segments = ([Segment([helix], exit='C'), ] +
                [Segment([helix], 'N', 'C')] * 3 +
                [Segment([helix], entry='N')])
    worms = grow(segments, SegmentSym('C2', lever=20),
                 executor=ProcessPoolExecutor, max_workers=2)
    assert 0.1411 < np.min(worms.scores) < 0.1412
    assert np.sum(worms.scores < 0.1412) == 4


@pytest.mark.skipif('not rcl.HAVE_PYROSETTA')
def test_grow_errors(curved_helix_pose):
    nsplice = SpliceSite(sele=[1, 2, 3, 4, 5, 6], polarity='N')
    csplice = SpliceSite(sele=[13, ], polarity='C')
    spliceable1 = Spliceable(body=curved_helix_pose, sites=[nsplice, csplice])
    spliceable2 = Spliceable(body=curved_helix_pose, sites=[nsplice, csplice])
    spliceables = [spliceable1]
    segments = ([Segment(spliceables, exit='C'), ] +
                [Segment(spliceables, 'N', 'C'), ] * 3 +
                [Segment(spliceables, entry='N'), ])
    checkc3 = SegmentSym('C2', from_seg=0, to_seg=-1)

    # make sure incorrect begin/end throws error
    with pytest.raises(ValueError):
        grow(segments[: 2], criteria=checkc3)
    with pytest.raises(ValueError):
        grow(segments[1:], criteria=checkc3)
    segments_polarity_mismatch = [
        Segment(spliceables, exit='C'),
        Segment(spliceables, entry='C'),
    ]
    with pytest.raises(ValueError):
        grow(segments_polarity_mismatch, criteria=checkc3)


@pytest.mark.skipif('not rcl.HAVE_PYROSETTA')
def test_memlimit(curved_helix_pose):
    helix = Spliceable(curved_helix_pose, sites=[((1, 2), 'N'), ('-2:', 'C')])
    segments = ([Segment([helix], exit='C'), ] +
                [Segment([helix], 'N', 'C')] * 3 +
                [Segment([helix], entry='N')])
    for i in range(2, 7):
        w1 = grow(segments, SegmentSym('c2'), memlim=10**i, thresh=30)
        assert i == 2 or np.allclose(w0.scores, w1.scores)
        w0 = w1


@pytest.mark.skipif('not rcl.HAVE_PYROSETTA')
def test_pose_alignment_0(curved_helix_pose):
    helix = Spliceable(curved_helix_pose, sites=[(1, 'N'), ('-4:', 'C')])
    segments = ([Segment([helix], exit='C'), ] +
                [Segment([helix], 'N', 'C')] * 3 +
                [Segment([helix], entry='N')])
    w = grow(segments, SegmentSym('c2'), thresh=1)
    assert len(w)
    assert tuple(w.indices[0]) in ((2, 1, 2, 0, 0), (1, 2, 0, 2, 0))
    pose = w.pose(0, align=1, withend=1)
    xyz0 = np.array([pose.residue(1).xyz(2)[i] for i in (0, 1, 2)] + [1])
    # resid 43 happens to be the symmetrically related one for this solution
    xyz1 = np.array([pose.residue(42).xyz(2)[i] for i in (0, 1, 2)] + [1])
    xyz1 = hrot([0, 0, 1], 180) @ xyz1
    assert np.sum((xyz1 - xyz0)**2) < 0.1


def show_with_axis(worms, idx=0):
    pose = worms.pose(idx, align=0, withend=1)
    showme(pose)
    import pymol
    pymol.finish_launching()
    x_from = worms.positions[idx][worms.criteria[0].from_seg]
    x_to = worms.positions[idx][worms.criteria[0].to_seg]
    x = x_to @ inv(x_from)
    axis, ang, cen = axis_ang_cen_of(x)
    print(ang)
    axis *= 100
    showline(axis, cen)
    showsphere(cen)


@pytest.mark.skipif('not rcl.HAVE_PYROSETTA')
def test_pose_alignment_1(curved_helix_pose):
    helix = Spliceable(curved_helix_pose, sites=[(1, 'N'), ('-4:', 'C')])
    segments = ([Segment([helix], exit='C'), ] +
                [Segment([helix], 'N', 'C')] * 4 +
                [Segment([helix], entry='N')])
    c = SegmentSym('c2', from_seg=1)
    w = grow(segments, c, thresh=1)
    assert len(w)
    # show_with_axis(w)
    # return
    pose = w.pose(0, align=1, withend=1)
    xyz0 = np.array([pose.residue(14).xyz(2)[i] for i in (0, 1, 2)] + [1])
    # resid 43 happens to be the symmetrically related one for this solution
    xyz1 = np.array([pose.residue(55).xyz(2)[i] for i in (0, 1, 2)] + [1])
    xyz1 = hrot([0, 0, 1], 180) @ xyz1
    assert np.sum((xyz1 - xyz0)**2) < 0.1


@pytest.mark.skipif('not rcl.HAVE_PYROSETTA')
def test_body_same_as_last(curved_helix_pose):
    helix = Spliceable(curved_helix_pose, sites=[(1, 'N'), ('-4:', 'C')])
    segments = ([Segment([helix, helix], exit='C'), ] +
                [Segment([helix], 'N', 'C')] * 3 +
                [Segment([helix, helix], entry='N')])
    w = grow(segments, SegmentSym('c2'), thresh=1)
    assert len(w)
    assert tuple(w.indices[0]) in ((2, 1, 2, 0, 0), (1, 2, 0, 2, 0))
    print(w.scores[0])
    for i, s in zip(w.indices, w.scores):
        assert segments[0].bodyid[i[0]] == segments[-1].bodyid[i[-1]]

    segments = ([Segment([helix], exit='C'), ] +
                [Segment([helix, helix], 'N', 'C'), ] +
                [Segment([helix], 'N', 'C')] * 3 +
                [Segment([helix, helix], entry='N')])
    w = grow(segments, SegmentSym('c2', 1), thresh=1)
    assert len(w)
    for i, s in zip(w.indices, w.scores):
        assert segments[1].bodyid[i[1]] == segments[-1].bodyid[i[-1]]
    assert tuple(w.indices[0]) in ((2, 2, 0, 2, 1, 0),
                                   (3, 2, 0, 2, 1, 0),
                                   (3, 6, 0, 2, 1, 1))


def test_reorder_spliced_as_N_to_C():
    Q = reorder_spliced_as_N_to_C

    with pytest.raises(ValueError): Q([[1], [1], [1]], 'NC')
    with pytest.raises(ValueError): Q([[1], [1], [1]], 'CN')
    with pytest.raises(ValueError): Q([[1, 1], [1], [1, 1]], 'CN')
    with pytest.raises(ValueError): Q([], 'CN')
    with pytest.raises(ValueError): Q([], '')
    with pytest.raises(ValueError): Q([[]], '')

    assert Q([[1]], '') == [[1]]
    assert Q([[1, 2]], '') == [[1], [2]]
    assert Q([[1], [2]], 'N') == [[1, 2]]
    assert Q([[1, 2], [3]], 'N') == [[1], [2, 3]]
    assert Q([[1, 2], [3, 4]], 'N') == [[1], [2, 3], [4]]
    assert Q([[1, 2, 3], [4, 5]], 'N') == [[1], [2], [3, 4], [5]]
    assert Q([[1], [2]], 'C') == [[2, 1]]
    assert Q([[1, 2], [3]], 'C') == [[1], [3, 2]]
    assert Q([[1, 2], [3, 4]], 'C') == [[1], [3, 2], [4]]
    assert Q([[1, 2, 3], [4, 5]], 'C') == [[1], [2], [4, 3], [5]]

    assert Q([[1], [2], [3]], 'NN') == [[1, 2, 3]]
    assert Q([[1], [2], [3, 4]], 'NN') == [[1, 2, 3], [4]]
    assert Q([[1], [2, 3], [4, 5]], 'NN') == [[1, 2], [3, 4], [5]]
    assert Q([[1, 2], [3, 4], [5, 6]], 'NN') == [[1], [2, 3], [4, 5], [6]]
    assert (Q([[1, 2, 3], [4, 5, 6], [7, 8, 9]], 'NN')
            == [[1], [2], [3, 4], [5], [6, 7], [8], [9]])
    assert (Q([[1, 2, 3], [4, 5, 6], [7, 8, 9]], 'CN')
            == [[1], [2], [4, 3], [5], [6, 7], [8], [9]])
    assert (Q([[1, 2, 3], [4, 5, 6], [7, 8, 9]], 'CC')
            == [[1], [2], [4, 3], [5], [7, 6], [8], [9]])
    assert (Q([[1, 2, 3], [4, 5, 6], [7, 8, 9]], 'NC')
            == [[1], [2], [3, 4], [5], [7, 6], [8], [9]])

    for n in range(10):
        x = [[i] for i in range(n + 1)]
        y = list(range(n + 1))
        assert Q(x, 'N' * n) == [y]
        assert Q(x, 'C' * n) == [y[::-1]]
        assert Q([[13, 14]] + x, 'N' + 'N' * n) == [[13], [14] + y]
        assert Q([[13, 14]] + x, 'C' + 'C' * n) == [[13], y[::-1] + [14]]
        assert (Q([[10, 11, 12]] + x + [[13, 14, 15]], 'N' + 'N' * n + 'N')
                == [[10], [11], [12] + y + [13], [14], [15]])
        assert (Q([[10, 11, 12]] + x + [[13, 14, 15]], 'C' + 'C' * n + 'C')
                == [[10], [11], [13] + y[::-1] + [12], [14], [15]])

    assert (Q([[1, 2, 3], [4, 5, 6], [7, 8, 9], [0, 1, 2]], 'NNN')
            == [[1], [2], [3, 4], [5], [6, 7], [8], [9, 0], [1], [2]])
    assert (Q([[1, 2, 3], [4, 5, 6], [7, 8, 9], [0, 1, 2]], 'CNN')
            == [[1], [2], [4, 3], [5], [6, 7], [8], [9, 0], [1], [2]])
    assert (Q([[1, 2, 3], [4, 5, 6], [7, 8, 9], [0, 1, 2]], 'NCN')
            == [[1], [2], [3, 4], [5], [7, 6], [8], [9, 0], [1], [2]])
    assert (Q([[1, 2, 3], [4, 5, 6], [7, 8, 9], [0, 1, 2]], 'NNC')
            == [[1], [2], [3, 4], [5], [6, 7], [8], [0, 9], [1], [2]])
    assert (Q([[1, 2, 3], [4, 5, 6], [7, 8, 9], [0, 1, 2]], 'NCC')
            == [[1], [2], [3, 4], [5], [7, 6], [8], [0, 9], [1], [2]])
    assert (Q([[1, 2, 3], [4, 5, 6], [11], [7, 8, 9], [0, 1, 2]], 'NCCC')
            == [[1], [2], [3, 4], [5], [7, 11, 6], [8], [0, 9], [1], [2]])
    assert (Q([[1, 2, 3], [4, 5, 6], [11], [12], [7, 8, 9], [0, 1, 2]], 'NCCCN')
            == [[1], [2], [3, 4], [5], [7, 12, 11, 6], [8], [9, 0], [1], [2]])
    assert (Q([[1, 2, 5, 5, 3], [4, 5, 6], [11], [12],
               [7, 8, 9], [0, 1, 2]], 'NCCCN')
            == [[1], [2], [5], [5], [3, 4], [5], [7, 12, 11, 6],
                [8], [9, 0], [1], [2]])


# @pytest.mark.xfail
@pytest.mark.skipif('not rcl.HAVE_PYROSETTA')
def test_make_pose_chains_dimer(c2pose):
    dimer = Spliceable(c2pose, sites=[('1,2:2', 'N'), ('2,3:3', 'N'),
                                      ('1,-4:-4', 'C'), ('2,-5:-5', 'C')])
    print(dimer)
    seq = dimer.body.sequence()[:12]

    dimerseg = Segment([dimer], 'N', '')
    enex, rest = dimerseg.make_pose_chains(0, pad=(0, 1))
    assert [x.sequence() for x in enex] == [seq[1:], seq]
    assert [x.sequence() for x in rest] == []
    assert enex[-1] is dimer.chains[2]
    enex, rest = dimerseg.make_pose_chains(1, pad=(0, 1))
    assert [x.sequence() for x in enex] == [seq[2:], seq]
    assert [x.sequence() for x in rest] == []
    assert enex[-1] is dimer.chains[1]

    dimerseg = Segment([dimer], 'C', '')
    enex, rest = dimerseg.make_pose_chains(0, pad=(0, 1))
    assert [x.sequence() for x in enex] == [seq[:-3], seq]
    assert [x.sequence() for x in rest] == []
    assert enex[-1] is dimer.chains[2]
    enex, rest = dimerseg.make_pose_chains(1, pad=(0, 1))
    assert [x.sequence() for x in enex] == [seq[:-4], seq]
    assert [x.sequence() for x in rest] == []
    assert enex[-1] is dimer.chains[1]

    dimerseg = Segment([dimer], '', 'N')
    enex, rest = dimerseg.make_pose_chains(0, pad=(0, 1))
    assert [x.sequence() for x in enex] == [seq, seq[1:]]
    assert [x.sequence() for x in rest] == []
    assert enex[0] is dimer.chains[2]
    enex, rest = dimerseg.make_pose_chains(1, pad=(0, 1))
    assert [x.sequence() for x in enex] == [seq, seq[2:]]
    assert [x.sequence() for x in rest] == []
    assert enex[0] is dimer.chains[1]

    dimerseg = Segment([dimer], 'N', 'N')
    enex, rest = dimerseg.make_pose_chains(0, pad=(0, 1))
    assert [x.sequence() for x in enex] == [seq[1:], seq[2:]]
    assert [x.sequence() for x in rest] == []
    enex, rest = dimerseg.make_pose_chains(1, pad=(0, 1))
    assert [x.sequence() for x in enex] == [seq[2:], seq[1:]]
    assert [x.sequence() for x in rest] == []
    with pytest.raises(IndexError):
        enex, rest = dimerseg.make_pose_chains(2, pad=(0, 1))

    dimerseg = Segment([dimer], 'N', 'C')
    enex, rest = dimerseg.make_pose_chains(0, pad=(0, 1))
    assert [x.sequence() for x in enex] == [seq[1:-3]]
    assert [x.sequence() for x in rest] == [seq]
    assert rest[0] is dimer.chains[2]
    enex, rest = dimerseg.make_pose_chains(1, pad=(0, 1))
    assert [x.sequence() for x in enex] == [seq[1:], seq[:-4]]
    assert [x.sequence() for x in rest] == []
    enex, rest = dimerseg.make_pose_chains(2, pad=(0, 1))
    assert [x.sequence() for x in enex] == [seq[2:], seq[:-3]]
    assert [x.sequence() for x in rest] == []
    enex, rest = dimerseg.make_pose_chains(3, pad=(0, 1))
    assert [x.sequence() for x in enex] == [seq[2:-4]]
    assert [x.sequence() for x in rest] == [seq]
    assert rest[0] is dimer.chains[1]
    with pytest.raises(IndexError):
        enex, rest = dimerseg.make_pose_chains(4, pad=(0, 1))


def residue_coords(p, ir, n=3):
    crd = (p.residue(ir).xyz(i) for i in range(1, n + 1))
    return np.stack([c.x, c.y, c.z, 1] for c in crd)


def residue_sym_err(p, ang, ir, jr, n=1):
    mxdist = 0
    for i in range(n):
        xyz0 = residue_coords(p, ir + i)
        xyz1 = residue_coords(p, jr + i)
        xyz3 = hrot([0, 0, 1], ang) @ xyz1.T
        xyz4 = hrot([0, 0, 1], -ang) @ xyz1.T
        mxdist = max(mxdist, min(
            np.max(np.sum((xyz0 - xyz3.T)**2, axis=1)),
            np.max(np.sum((xyz0 - xyz4.T)**2, axis=1))))
    return np.sqrt(mxdist)


@pytest.mark.skipif('not rcl.HAVE_PYROSETTA')
def test_multichain_reveres_pol(c2pose, curved_helix_pose):
    helix = Spliceable(
        curved_helix_pose, sites=[((1, 3,), 'N'), ((10, 11, 13), 'C')])
    dimer = Spliceable(c2pose, sites=[('1,:2', 'N'), ('1,-1:', 'C'),
                                      ('2,:2', 'N'), ('2,-1:', 'C')])
    segments = [Segment([helix], exit='C'),
                Segment([helix], entry='N', exit='C'),
                Segment([dimer], entry='N', exit='C'),
                Segment([helix], entry='N', exit='C'),
                Segment([helix], entry='N'), ]
    wnc = grow(segments, SegmentSym('C3'), thresh=1)
    assert len(wnc)
    assert residue_sym_err(wnc.pose(0), 120, 2, 54, 8) < 0.5

    segments = [Segment([helix], exit='N'),
                Segment([helix], entry='C', exit='N'),
                Segment([dimer], entry='C', exit='N'),
                Segment([helix], entry='C', exit='N'),
                Segment([helix], entry='C'), ]
    wcn = grow(segments, SegmentSym('C3'), thresh=1)
    assert residue_sym_err(wcn.pose(0), 120, 22, 35, 8) < 0.5

    # N-to-C and C-to-N construction should be same
    assert np.allclose(wnc.scores, wcn.scores, atol=1e-3)


@pytest.mark.skipif('not rcl.HAVE_PYROSETTA')
def test_multichain_mixed_pol(c2pose, c3pose, curved_helix_pose):
    helix = Spliceable(curved_helix_pose, [(':4', 'N'), ((10, 12, 13), 'C')])
    dimer = Spliceable(c2pose, sites=[('1,:2', 'N'), ('1,-1:', 'C'),
                                      ('2,:2', 'N'), ('2,-1:', 'C')])
    trimer = Spliceable(c3pose, sites=[('1,:1', 'N'), ('1,-2:', 'C'),
                                       # ('2,:2', 'N'), ('2,-2:', 'C'),
                                       ('3,:1', 'N'), ('3,-2:', 'C')])
    segments = [Segment([helix], exit='C'),
                Segment([dimer], entry='N', exit='N'),
                Segment([helix], entry='C', exit='N'),
                Segment([trimer], entry='C', exit='C'),
                Segment([helix], entry='N')]
    w = grow(segments, SegmentSym('C3'), thresh=2)
    assert residue_sym_err(w.pose(0), 120, 2, 67, 9) < 2
    w = grow(segments, SegmentSym('C4'), thresh=2)
    assert residue_sym_err(w.pose(0), 90, 2, 69, 9) < 2
    w = grow(segments, SegmentSym('C5'), thresh=2)
    assert residue_sym_err(w.pose(0), 72, 2, 65, 9) < 2
    w = grow(segments, SegmentSym('C6'), thresh=2)
    assert residue_sym_err(w.pose(0), 60, 2, 71, 9) < 2


@pytest.mark.skipif('not rcl.HAVE_PYROSETTA')
def test_multichain_db(c2pose, c3pose, curved_helix_pose):
    helix = Spliceable(curved_helix_pose, [(':4', 'N'), ('-4:', "C")])
    dimer = Spliceable(c2pose, sites=[('1,-1:', 'C'), ('1,-1:', 'C')])
    segments = [Segment([helix], exit='N'),
                Segment([dimer], entry='C', exit='C'),
                Segment([helix], entry='N')]
    w = grow(segments, SegmentSym('C4'), thresh=20)
    with pytest.raises(ValueError):
        showme(w.pose(0))
