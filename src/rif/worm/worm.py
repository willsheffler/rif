from rif.rcl import Pose, bbstubs, to_rosetta_stub
import rosetta
from rif.homo import axis_angle_of, angle_of, hrot
from rif.vis import showme
import numpy as np
from numpy.linalg import inv, norm
import math

identity44f4 = np.identity(4, dtype='f4')
identity44f8 = np.identity(4, dtype='f8')


class WormCriteria:
    pass


class AxesIntersect(WormCriteria):
    pass


class SegmentXform(WormCriteria):

    def __init__(self, symmetry, *, tol=1.0, from_seg=0,
                 orgnin_seg=None, lever=100.0, to_seg=-1):
        self.symmetry = symmetry
        self.tol = tol
        self.from_seg = from_seg
        self.orgnin_seg = orgnin_seg
        self.lever = lever
        self.to_seg = to_seg
        self.rot_tol = tol / lever
        if self.symmetry[0] in 'cC':
            self.nfold = int(self.symmetry[1:])
            if self.nfold <= 0:
                raise ValueError('invalid symmetry: ' + symmetry)
            self.symangle = math.pi * 2.0 / self.nfold
        else: raise ValueError('can only do Cx symmetry for now')
        assert not orgnin_seg
        if self.tol <= 0: raise ValueError('tol should be > 0')

    def score(self, segment_pos, *, debug=False, **kwarge):
        x_from = segment_pos[self.from_seg]
        x_to = segment_pos[self.to_seg]
        xhat = inv(x_from) @ x_to
        trans = xhat[..., :3, 3]
        if self.orgnin_seg:
            raise NotImplementedError
        elif self.nfold is 1:
            angle = angle_of(xhat)
            carterrsq = np.sum(trans**2, axis=-1)
            roterrsq = angle**2
        else:
            axis, angle = axis_angle_of(xhat)
            carterrsq = np.sum(trans * axis, axis=-1)**2
            roterrsq = (angle - self.symangle)**2
        return np.sqrt(carterrsq / self.tol**2 + roterrsq / self.rot_tol**2)


class SpliceSite:

    def __init__(self, sele, polarity):
        if isinstance(sele, str) or isinstance(sele, int):
            sele = [sele]
        self.selections = list(sele)
        self.polarity = polarity

    def resid(self, id, body):
        resid = id if id >= 0 else len(body) + 1 + id
        if not 0 < resid <= len(body):
            raise ValueError('resid ' + str(resid)
                             + ' invalid for body of size ' + str(len(body)))
        return resid

    def resids(self, body):
        resids = set()
        for sele in self.selections:
            try:
                if isinstance(sele, int):
                    resids.add(self.resid(sele, body))
                elif isinstance(sele, str):
                    s = sele.split(':')
                    start = self.resid(int(s[0] or 1), body)
                    stop = self.resid(int(s[1] or -1), body)
                    step = int(s[2]) if len(s) > 2 else 1
                    # print(start, stop + 1, step)
                    for ir in range(start, stop + 1, step):
                        resids.add(ir)
                elif sele is None:
                    resids.add(None)
                else:
                    raise ValueError('selection must be int, str, or None')
            except ValueError as e:
                raise ValueError('Error with selection '
                                 + str(sele) + ': ' + str(e))
        resids = sorted(resids)
        if not resids:
            raise ValueError('empty SpliceSite')
        return resids


class Spliceable:

    def __init__(self, body, *, sites, bodyid=None):
        self.body = body
        self.bodyid = bodyid
        if callable(sites):
            sites = sites(body)
        if isinstance(sites, SpliceSite):
            sites = [sites]
        self.sites = list(sites)
        for i in range(len(self.sites)):
            if not isinstance(self.sites[i], SpliceSite):
                assert len(self.sites[i]) is 2
                self.sites[i] = SpliceSite(*self.sites[i])

    def splicable_positions(self):
        """selection of resids, and map 'global' index to selected index"""
        resid_subset = set()
        for site in self.sites:
            resid_subset |= set(site.resids(self.body))
        resid_subset = np.array(list(resid_subset))
        # really? must be an easier way to 'invert' a mapping in numpy?
        N = len(self.body) + 1
        val, idx = np.where(0 == (np.arange(N)[np.newaxis, :] -
                                  resid_subset[:, np.newaxis]))
        to_subset = np.array(N * [-1])
        to_subset[idx] = val
        assert np.all(to_subset[resid_subset] == np.arange(len(resid_subset)))
        return resid_subset, to_subset


class Segment:

    def __init__(self, splicables, *, entry=None, exit=None):
        self.entrypol = entry
        self.exitpol = exit
        self.init(splicables, entry, exit)

    def init(self, splicables=None, entry=None, exit=None):
        if not (entry or exit):
            raise ValueError('at least one of entry/exit required')
        self.splicables = list(splicables) or self.splicables
        self.entrypol = entry or self.entrypol
        self.exitpol = exit or self.exitpol
        # each array has all in/out pairs
        self.x2exit, self.x2orgn = list(), list()
        self.entryresid, self.exitresid, self.bodyid = list(), list(), list()
        # this whole loop is pretty inefficient, but that probably
        # doesn't matter much given the cost subsequent operations (?)
        for ibody, splicable in enumerate(self.splicables):
            resid_subset, to_subset = splicable.splicable_positions()
            bodyid = ibody if splicable.bodyid is None else splicable.bodyid
            # extract 'stubs' from body at selected positions
            # rif 'stubs' have 'extra' 'features'... the raw field is
            # just bog-standard homogeneous matrices
            stubs = bbstubs(splicable.body, resid_subset)['raw']
            if len(resid_subset) != stubs.shape[0]:
                raise ValueError("no funny residues supported")
            stubs_inv = inv(stubs)
            entry_sites = (list(enumerate(splicable.sites)) if self.entrypol else
                           [(-1, SpliceSite(sele=[None], polarity=self.entrypol))])
            exit_sites = (list(enumerate(splicable.sites)) if self.exitpol else
                          [(-1, SpliceSite(sele=[None], polarity=self.exitpol))])
            for isite, entry_site in entry_sites:
                if entry_site.polarity == self.entrypol:
                    for jsite, exit_site in exit_sites:
                        if isite != jsite and exit_site.polarity == self.exitpol:
                            for ires in entry_site.resids(splicable.body):
                                istub_inv = (identity44f4 if not ires
                                             else stubs_inv[to_subset[ires]])
                                ires = ires or -1
                                for jres in exit_site.resids(splicable.body):
                                    jstub = (identity44f4 if not jres
                                             else stubs[to_subset[jres]])
                                    jres = jres or -1
                                    self.x2exit.append(istub_inv @ jstub)
                                    self.x2orgn.append(istub_inv)
                                    self.entryresid.append(ires)
                                    self.exitresid.append(jres)
                                    self.bodyid.append(bodyid)
        if len(self.x2exit) is 0:
            raise ValueError('no valid splices found')
        self.x2exit = np.stack(self.x2exit)
        self.x2orgn = np.stack(self.x2orgn)
        self.entryresid = np.array(self.entryresid)
        self.exitresid = np.array(self.exitresid)
        self.bodyid = np.array(self.bodyid)

    def make_pose(self, index, append_to=None, position=None):
        append_to = append_to or Pose()
        pose = self.splicables[self.bodyid[index]].body
        lb = self.entryresid[index]
        ub = self.exitresid[index]
        # print('make_pose', lb, ub)
        lb = lb if lb > 0 else 1
        ub = ub if ub > 0 else len(pose)
        if append_to and self.exitresid[index] > 0: ub = ub - 1
        # lb, ub = 1, len(pose) + 1

        assert 0 < lb <= ub - 1 <= len(pose)
        # print('append_subpose', lb, ub - 1)
        rosetta.core.pose.append_subpose_to_pose(append_to, pose, lb, ub)

        if position is not None:
            x = to_rosetta_stub(position)
            xub = len(append_to)
            xlb = xub - ub + lb
            assert 0 < xlb <= xub <= len(append_to)
            # print('xform_pose', xlb, xub)
            rosetta.protocols.sic_dock.xform_pose(append_to, x, xlb, xub)
        return append_to


class Worms:

    def __init__(self, segments, score, solutions):
        self.segments = segments
        self.score = score
        self.solutions = solutions


def chain_xforms(x2exit, x2orgn):
    fullaxes = (np.newaxis,) * (len(x2exit) - 1)
    xconn = [x2exit[0][fullaxes], ]
    xbody = [x2orgn[0][fullaxes], ]
    for iseg in range(1, len(x2exit)):
        fullaxes = (slice(None),) + (np.newaxis,) * iseg
        xconn.append(xconn[iseg - 1] @ x2exit[iseg][fullaxes])
        xbody.append(xconn[iseg - 1] @ x2orgn[iseg][fullaxes]
                     if iseg + 1 < len(x2exit) else xconn[iseg])
    perm = list(range(len(xbody) - 1, -1, -1)) + [len(xbody), len(xbody) + 1]
    xbody = [np.transpose(x, perm) for x in xbody]
    xconn = [np.transpose(x, perm) for x in xconn]
    return xbody, xconn


def grow(segments, criteria, *, cache=None, threshold=1):
    if segments[0].entrypol is not None:
        raise ValueError('beginning of worm cant have entry')
    if segments[-1].exitpol is not None:
        raise ValueError('end of worm cant have exit')
    criteria = [criteria] if isinstance(criteria, WormCriteria) else criteria
    for a, b in zip(segments[:-1], segments[1:]):
        if not (a.exitpol and b.entrypol and a.exitpol != b.entrypol):
            raise ValueError('incompatible exit->entry polarity: '
                             + str(a.exitpol) + '->'
                             + str(b.entrypol) + ' on segment pair: '
                             + str((segments.index(a), segments.index(b))))

    sizes = [len(s.bodyid) for s in segments]
    print('grow: size', sizes, np.prod(sizes))

    segment_pos, connect_pos = chain_xforms([s.x2exit for s in segments],
                                            [s.x2orgn for s in segments])
    assert segment_pos[-1].shape[:-2] == tuple(len(s.bodyid) for s in segments)

    score = sum(c.score(segment_pos=segment_pos) for c in criteria)
    imin = np.unravel_index(np.argmin(score), score.shape)

    print('best:', score[imin], imin)

    # return

    poses = list()
    # poses = Pose()
    positions = list()
    for iseg, seg in enumerate(segments):
        i = imin[iseg]
        pos = segment_pos[iseg]
        # print('segment %2i body %3i entry %-04i exit %-04i' % (
        # iseg, seg.bodyid[i], seg.entryresid[i], seg.exitresid[i]))
        ipos = tuple(np.minimum(np.array(pos.shape[:-2]) - 1, imin))
        position = pos[ipos]
        positions.append(position)
        poses.append(seg.make_pose(i, position=position))
        # poses = seg.make_pose(i, position=position, append_to=poses)

    score2 = criteria[0].score(positions, debug=1)
    print('rescore', score2)

    return

    print(poses)
    showme(poses)
    from pymol import cmd
    cmd.hide('ev')
    cmd.show('lines', 'name n+ca+c')
    import time
    while 1:
        time.sleep(1)

    return None
