from rif.worm import *
from rif.vis import showme
from rif import rcl
import pytest
import numpy as np


def test_geom_check():
    with pytest.raises(ValueError):
        GeomCheck()
    with pytest.raises(ValueError):
        GeomCheck(is_exactly='dummy', sym='dummy')
    with pytest.raises(ValueError):
        GeomCheck(is_exactly='dummy', z_axes_intersection='dummy')
    with pytest.raises(ValueError):
        GeomCheck(is_exactly='dummy', origin_segment='dummy')
    with pytest.raises(ValueError):
        GeomCheck(is_exactly='dummy', origin_segment='dummy')


@pytest.mark.skipif('not rcl.HAVE_PYROSETTA')
def test_segment_geom(curved_helix_pose):
    "currently only a basic sanity checkb... only checks translation distances"
    body = curved_helix_pose
    nsplice = SpliceSite(polarity='N', resids=[1, 2, ])
    csplice = SpliceSite(polarity='C', resids=[9, 10, 11, 12, 13])
    Npairs0 = len(nsplice.resids) * len(csplice.resids)

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
    assert seg.x2exit.shape == (Nexsite * len(csplice.resids), 4, 4)
    assert seg.x2orig.shape == (Nexsite * len(csplice.resids), 4, 4)
    assert np.all(seg.x2exit[..., 3, :3] == 0)
    assert np.all(seg.x2exit[..., 3, 3] == 1)
    for e2x, e2o, ir, jr in zip(seg.x2exit, seg.x2orig,
                                seg.entryresid, seg.exitresid):
        assert np.isnan(ir)
        jdis = np.sqrt(np.sum((e2x[..., :3, 3]**2)))
        cadis = body.residue(jr).xyz('CA').length()
        assert np.all(e2o == np.identity(4))
        assert abs(jdis - cadis) < 0.001

    # test middle segment with entry and exit
    seg = Segment([splicable], entry='N', exit='C')
    assert seg.x2exit.shape == (Nexsite**2 * Npairs0, 4, 4)
    assert seg.x2orig.shape == (Nexsite**2 * Npairs0, 4, 4)
    assert np.all(seg.x2exit[..., 3, :3] == 0)
    assert np.all(seg.x2exit[..., 3, 3] == 1)
    for e2x, e2o, ir, jr in zip(seg.x2exit, seg.x2orig,
                                seg.entryresid, seg.exitresid):
        ijdis = np.sqrt(np.sum((e2x[..., :3, 3]**2)))
        cadis = body.residue(ir).xyz('CA').distance(body.residue(jr).xyz('CA'))
        assert abs(ijdis - cadis) < 0.001
        odis = np.sqrt(np.sum((e2o[..., :3, 3]**2)))
        caodis = body.residue(ir).xyz('CA').length()
        assert abs(odis - caodis) < 0.001

    # test ending segment.. only has entry
    seg = Segment([splicable], entry='N')
    assert seg.x2exit.shape == (Nexsite * len(nsplice.resids), 4, 4)
    assert seg.x2orig.shape == (Nexsite * len(nsplice.resids), 4, 4)
    assert np.all(seg.x2exit[..., 3, :3] == 0)
    assert np.all(seg.x2exit[..., 3, 3] == 1)
    for e2x, e2o, ir, jr in zip(seg.x2exit, seg.x2orig,
                                seg.entryresid, seg.exitresid):
        assert np.isnan(jr)
        idis = np.sqrt(np.sum((e2x[..., :3, 3]**2)))
        cadis = body.residue(ir).xyz('CA').length()
        assert abs(idis - cadis) < 0.001
        assert np.all(e2o == e2x)

    # test now with multiple splicables input to segment
    Nexbody = 3
    seg = Segment([splicable] * Nexbody, entry='N', exit='C')
    Npairs_expected = Nexbody * Nexsite**2 * Npairs0
    assert seg.x2exit.shape == (Npairs_expected, 4, 4)
    assert seg.x2orig.shape == (Npairs_expected, 4, 4)
    assert len(seg.entryresid) == Npairs_expected
    assert len(seg.exitresid) == Npairs_expected
    assert len(seg.bodyid) == Npairs_expected
    for i in range(Nexbody):
        assert i == seg.bodyid[0 + i * Npairs0 * Nexsite**2]
    assert np.all(seg.x2exit[..., 3, :3] == 0)
    assert np.all(seg.x2exit[..., 3, 3] == 1)
    for e2x, e2o, ir, jr in zip(seg.x2exit, seg.x2orig,
                                seg.entryresid, seg.exitresid):
        ijdis = np.sqrt(np.sum((e2x[..., :3, 3]**2)))
        cadis = body.residue(ir).xyz('CA').distance(body.residue(jr).xyz('CA'))
        assert abs(ijdis - cadis) < 0.001
        odis = np.sqrt(np.sum((e2o[..., :3, 3]**2)))
        caodis = body.residue(ir).xyz('CA').length()
        assert abs(odis - caodis) < 0.001


# todo: move elsewhere???
@pytest.mark.skipif('not rcl.HAVE_PYROSETTA')
def test_all_xform_combinations(curved_helix_pose):
    pass


@pytest.mark.skipif('not rcl.HAVE_PYROSETTA')
def test_exmple(curved_helix_pose):
    nsplice = SpliceSite(resids=[1, 2], polarity='N')
    csplice = SpliceSite(resids=[11, 12, 13, ], polarity='C')
    splicable = Spliceable(body=curved_helix_pose, sites=[nsplice, csplice])
    splicables = [splicable, ]
    segments = [Segment(splicables, exit='C'),
                Segment(splicables, entry='N', exit='C'),
                Segment(splicables, entry='N')]
    checkc3 = GeomCheck(from_segment=0, to_segment=-1, sym='C3')
    grow(segments, criteria=checkc3)


@pytest.mark.skipif('not rcl.HAVE_PYROSETTA')
def test_grow(curved_helix_pose, N=6):
    nsplice = SpliceSite(resids=[1, 2, ], polarity='N')
    csplice = SpliceSite(resids=[9, 10, 11, 12, 13, ], polarity='C')
    splicable1 = Spliceable(body=curved_helix_pose, sites=[nsplice, csplice])
    splicable2 = Spliceable(body=curved_helix_pose, sites=[nsplice, csplice])
    splicables = [splicable1]
    for i in range(N):
        segments = ([Segment(splicables, exit='C'), ] +
                    [Segment(splicables, entry='N', exit='C'), ] * i +
                    [Segment(splicables, entry='N'), ])
        checkc3 = GeomCheck(from_segment=0, to_segment=-1, sym='C3')
        grow(segments, criteria=checkc3)

    assert 0  # temporary!!!

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


if __name__ == '__main__':
    import pyrosetta
    pyrosetta.init()
    pose = pyrosetta.pose_from_file(
        '/home/sheffler/rifsrc/src/rif/data/pdb/curved_helix.pdb')
    print('foo')
    test_grow(pose, 30)
