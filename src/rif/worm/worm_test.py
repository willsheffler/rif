from rif import worm
from rif.vis import showme
from rif import rcl
import pytest
import numpy as np


def test_geom_check():
    with pytest.raises(ValueError):
        worm.GeomCheck()
    with pytest.raises(ValueError):
        worm.GeomCheck(is_exactly='dummy', sym='dummy')
    with pytest.raises(ValueError):
        worm.GeomCheck(is_exactly='dummy', z_axes_intersection='dummy')
    with pytest.raises(ValueError):
        worm.GeomCheck(is_exactly='dummy', origin_segment='dummy')
    with pytest.raises(ValueError):
        worm.GeomCheck(is_exactly='dummy', origin_segment='dummy')


@pytest.mark.skipif('not rcl.HAVE_PYROSETTA')
def test_segment_geom(curved_helix_pose):
    "currently only a basic sanity checkb... only checks translation distances"
    body = curved_helix_pose
    nsplice = worm.SpliceSite(polarity='N', resids=[1, 2, ])
    csplice = worm.SpliceSite(polarity='C', resids=[9, 10, 11, 12, 13])
    Npairs0 = len(nsplice.resids) * len(csplice.resids)

    # N to N and C to C invalid, can't splice to same
    splicable = worm.Splicable(body, [nsplice, csplice])
    with pytest.raises(ValueError):
        s = worm.Segment([splicable], entry='N', egress='N')
    with pytest.raises(ValueError):
        s = worm.Segment([splicable] * 3, entry='C', egress='C')

    # add some extra splice sites
    Nexsite = 2
    splicable = worm.Splicable(body, [nsplice, csplice] * Nexsite)

    # test beginning segment.. only has egress
    s = worm.Segment([splicable], egress='C')
    assert s.entry2exit.shape == (Nexsite * len(csplice.resids), 4, 4)
    assert s.entry2orig.shape == (Nexsite * len(csplice.resids), 4, 4)
    assert np.all(s.entry2exit[..., 3, :3] == 0)
    assert np.all(s.entry2exit[..., 3, 3] == 1)
    for e2x, e2o, ir, jr in zip(s.entry2exit, s.entry2orig,
                                s.entryresid, s.exitresid):
        assert np.isnan(ir)
        jdis = np.sqrt(np.sum((e2x[..., :3, 3]**2)))
        cadis = body.residue(jr).xyz('CA').length()
        assert np.all(e2o == np.identity(4))
        assert abs(jdis - cadis) < 0.001

    # test middle segment with entry and egress
    s = worm.Segment([splicable], entry='N', egress='C')
    assert s.entry2exit.shape == (Nexsite**2 * Npairs0, 4, 4)
    assert s.entry2orig.shape == (Nexsite**2 * Npairs0, 4, 4)
    assert np.all(s.entry2exit[..., 3, :3] == 0)
    assert np.all(s.entry2exit[..., 3, 3] == 1)
    for e2x, e2o, ir, jr in zip(s.entry2exit, s.entry2orig,
                                s.entryresid, s.exitresid):
        ijdis = np.sqrt(np.sum((e2x[..., :3, 3]**2)))
        cadis = body.residue(ir).xyz('CA').distance(body.residue(jr).xyz('CA'))
        assert abs(ijdis - cadis) < 0.001
        odis = np.sqrt(np.sum((e2o[..., :3, 3]**2)))
        caodis = body.residue(ir).xyz('CA').length()
        assert abs(odis - caodis) < 0.001

    # test ending segment.. only has entry
    s = worm.Segment([splicable], entry='N')
    assert s.entry2exit.shape == (Nexsite * len(nsplice.resids), 4, 4)
    assert s.entry2orig.shape == (Nexsite * len(nsplice.resids), 4, 4)
    assert np.all(s.entry2exit[..., 3, :3] == 0)
    assert np.all(s.entry2exit[..., 3, 3] == 1)
    for e2x, e2o, ir, jr in zip(s.entry2exit, s.entry2orig,
                                s.entryresid, s.exitresid):
        assert np.isnan(jr)
        idis = np.sqrt(np.sum((e2x[..., :3, 3]**2)))
        cadis = body.residue(ir).xyz('CA').length()
        assert abs(idis - cadis) < 0.001
        assert np.all(e2o == e2x)

    # test now with multiple splicables input to segment
    Nexbody = 3
    s = worm.Segment([splicable] * Nexbody, entry='N', egress='C')
    Npairs_expected = Nexbody * Nexsite**2 * Npairs0
    assert s.entry2exit.shape == (Npairs_expected, 4, 4)
    assert s.entry2orig.shape == (Npairs_expected, 4, 4)
    assert len(s.entryresid) == Npairs_expected
    assert len(s.exitresid) == Npairs_expected
    assert len(s.bodyid) == Npairs_expected
    for i in range(Nexbody):
        assert i == s.bodyid[0 + i * Npairs0 * Nexsite**2]
    assert np.all(s.entry2exit[..., 3, :3] == 0)
    assert np.all(s.entry2exit[..., 3, 3] == 1)
    for e2x, e2o, ir, jr in zip(s.entry2exit, s.entry2orig,
                                s.entryresid, s.exitresid):
        ijdis = np.sqrt(np.sum((e2x[..., :3, 3]**2)))
        cadis = body.residue(ir).xyz('CA').distance(body.residue(jr).xyz('CA'))
        assert abs(ijdis - cadis) < 0.001
        odis = np.sqrt(np.sum((e2o[..., :3, 3]**2)))
        caodis = body.residue(ir).xyz('CA').length()
        assert abs(odis - caodis) < 0.001


@pytest.mark.skipif('not rcl.HAVE_PYROSETTA')
def test_grow(curved_helix_pose):
    nsplice = worm.SpliceSite(resids=[1, 2, 3, ], polarity='N')
    csplice = worm.SpliceSite(resids=[10, 11, 12, 13], polarity='C')
    splicable = worm.Splicable(curved_helix_pose, [nsplice, csplice])
    segments = [
        worm.Segment([splicable], egress='C'),
        worm.Segment([splicable], entry='N', egress='C'),
        worm.Segment([splicable], entry='N'),
    ]
    checkc3 = worm.GeomCheck(from_segment=0, to_segment=-1, sym='C3')
    worm.grow(segments, criteria=checkc3)

    # make sure incorrect begin/end throws error
    with pytest.raises(ValueError):
        worm.grow(segments[:2], criteria=checkc3)
    with pytest.raises(ValueError):
        worm.grow(segments[1:], criteria=checkc3)
    segments_bad = [
        worm.Segment([splicable], egress='C'),
        worm.Segment([splicable], entry='C', egress='N'),
        worm.Segment([splicable], entry='N'),
    ]
    with pytest.raises(ValueError):
        worm.grow(segments_bad, criteria=checkc3)
