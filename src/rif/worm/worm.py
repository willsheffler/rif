from numpy.linalg import inv
from pyrosetta import rosetta as ros
from rif.homog import axis_angle_of, angle_of, hrot
from rif.rcl import Pose, bbstubs, to_rosetta_stub
from rif.vis import showme
from tqdm import tqdm
import functools as ft
import itertools as it
from concurrent.futures import ThreadPoolExecutor, as_completed
import multiprocessing
import numpy as np
import os
import sys

identity44f4 = np.identity(4, dtype='f4')
identity44f8 = np.identity(4, dtype='f8')


# todo: the following should go elsewhere...

def cpu_count():
    try:
        return int(os.environ['SLURM_CPUS_ON_NODE'])
    except:
        return multiprocessing.cpu_count()


def tqdm_parallel_map(pool, function, *args, **kw):
    futures = [pool.submit(function, *a) for a in zip(*args)]
    return (f.result() for f in tqdm(as_completed(futures),
                                     total=len(futures), **kw))


class WormCriteria: pass


class AxesIntersect(WormCriteria): pass


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
            self.symangle = np.pi * 2.0 / self.nfold
        else: raise ValueError('can only do Cx symmetry for now')
        assert not orgnin_seg
        if self.tol <= 0: raise ValueError('tol should be > 0')

    def score(self, segpos, *, debug=False, **kwarge):
        x_from = segpos[self.from_seg]
        x_to = segpos[self.to_seg]
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

    def __len__(self):
        return len(self.bodyid)

    def init(self, splicables=None, entry=None, exit=None):
        if not (entry or exit):
            raise ValueError('at least one of entry/exit required')
        self.splicables = list(splicables) or self.splicables
        self.entrypol = entry or self.entrypol
        self.exitpol = exit or self.exitpol
        # each array has all in/out pairs
        self.x2exit, self.x2orgn, self.bodyid = [], [], []
        self.entryresid, self.exitresid = [], []
        self.entrysiteid, self.exitsiteid = [], []
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
                                    self.entrysiteid.append(isite)
                                    self.entryresid.append(ires)
                                    self.exitsiteid.append(jsite)
                                    self.exitresid.append(jres)
                                    self.bodyid.append(bodyid)
        if len(self.x2exit) is 0:
            raise ValueError('no valid splices found')
        self.x2exit = np.stack(self.x2exit)
        self.x2orgn = np.stack(self.x2orgn)
        self.entrysiteid = np.stack(self.entrysiteid)
        self.entryresid = np.array(self.entryresid)
        self.exitsiteid = np.array(self.exitsiteid)
        self.exitresid = np.array(self.exitresid)
        self.bodyid = np.array(self.bodyid)

    def make_pose(self, index, append_to=None,
                  position=None, onechain=True):
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
        ros.core.pose.append_subpose_to_pose(
            append_to, pose, lb, ub, new_chain=not onechain)

        if position is not None:
            x = to_rosetta_stub(position)
            xub = len(append_to)
            xlb = xub - ub + lb
            assert 0 < xlb <= xub <= len(append_to)
            # print('xform_pose', xlb, xub)
            ros.protocols.sic_dock.xform_pose(append_to, x, xlb, xub)
        return append_to


class Worms:

    def __init__(self, segments, scores, indices, positions):
        self.segments = segments
        self.scores = scores
        self.indices = indices
        self.positions = positions

    def __init___(self):
        return len(self.scores)

    def pose(self, which, **kw):
        if hasattr(which, '__iter__'):
            return (self.pose(w) for w in which)
        pose = Pose()
        for iseg, seg in enumerate(self.segments):
            i = self.indices[which][iseg]
            position = self.positions[which][iseg]
            pose = seg.make_pose(i, position=position, append_to=pose, **kw)
        return pose

    def __len__(self):
        return len(self.scores)


def _chain_xforms(segments):
    x2exit = [s.x2exit for s in segments]
    x2orgn = [s.x2orgn for s in segments]
    fullaxes = (np.newaxis,) * (len(x2exit) - 1)
    xconn = [x2exit[0][fullaxes], ]
    xbody = [x2orgn[0][fullaxes], ]
    for iseg in range(1, len(x2exit)):
        fullaxes = (slice(None),) + (np.newaxis,) * iseg
        xconn.append(xconn[iseg - 1] @ x2exit[iseg][fullaxes])
        xbody.append(xconn[iseg - 1] @ x2orgn[iseg][fullaxes])
    perm = list(range(len(xbody) - 1, -1, -1)) + [len(xbody), len(xbody) + 1]
    xbody = [np.transpose(x, perm) for x in xbody]
    xconn = [np.transpose(x, perm) for x in xconn]
    return xbody, xconn


def _grow_chunk(samp, segpos, connpos, segs, end, criteria, thresh):
    segpos, connpos = segpos[:end], connpos[:end]
    for iseg, seg in enumerate(segs[end:]):
        segpos.append(connpos[-1] @ seg.x2orgn[samp[iseg]])
        connpos.append(connpos[-1] @ seg.x2exit[samp[iseg]])
    score = sum(c.score(segpos=segpos) for c in criteria)
    ilow0 = np.where(score < thresh)
    sampidx = tuple(np.repeat(i, len(ilow0[0])) for i in samp)
    lowpostmp = []
    for iseg in range(len(segpos)):
        ilow = ilow0[:iseg + 1] + (0,) * (segpos[0].ndim - 2 - (iseg + 1))
        lowpostmp.append(segpos[iseg][ilow])
    return score[ilow0], np.array(ilow0 + sampidx).T, np.stack(lowpostmp, 1)


def _grow_chunks(ijob, context):
    sampsizes, njob, segments, end, criteria, thresh = context
    samples = it.product(*(range(n) for n in sampsizes))
    segpos, connpos = _chain_xforms(segments[:end])  # common data
    args = [list(samples)[ijob::njob]] + [it.repeat(x) for x in (
        segpos, connpos, segments, end, criteria, thresh)]
    chunk = list(map(_grow_chunk, *args))
    return [np.concatenate([c[i] for c in chunk])
            for i in range(3)] if chunk else None


def grow(segments, criteria, *, last_body_same_as=None,
         cache=None, thresh=2, expert=False, memlim=1e6,
         executor=None, max_workers=None, debug=False, jobmult=32):
    # checks
    if segments[0].entrypol is not None:
        raise ValueError('beginning of worm cant have entry')
    if segments[-1].exitpol is not None:
        raise ValueError('end of worm cant have exit')
    if last_body_same_as is not None and not expert and (
            segments[last_body_same_as] is not segments[-1]):
        raise ValueError("segments[last_body_same_as] not same as segments[-1], "
                         + "if you're sure, pass expert=True")
    for a, b in zip(segments[:-1], segments[1:]):
        if not (a.exitpol and b.entrypol and a.exitpol != b.entrypol):
            raise ValueError('incompatible exit->entry polarity: '
                             + str(a.exitpol) + '->'
                             + str(b.entrypol) + ' on segment pair: '
                             + str((segments.index(a), segments.index(b))))
    if last_body_same_as is not None:
        raise NotImplementedError('must implement subselection for last seg')

    # setup
    criteria = [criteria] if isinstance(criteria, WormCriteria) else criteria
    if executor is None:
        executor = ThreadPoolExecutor
        max_workers = 1
    sizes = [len(s.bodyid) for s in segments]
    end = len(segments) - 1
    while end > 1 and (np.prod(sizes[end:]) < cpu_count() or
                       memlim <= 64 * np.prod(sizes[:end])): end -= 1
    ntot, nchunk, nchunks = (np.product(x)
                             for x in (sizes, sizes[:end], sizes[end:]))
    nworker = max_workers or cpu_count()
    njob = nworker * jobmult
    njob = min(njob, nchunks)
    print('tot = {:,}, nchunk = {:,}, nchunks = {:,}, nworker = {}, njob = {}, '
          'worm/job = {:,}, chunk/job = {}'.format(
              ntot, nchunk, nchunks, nworker, njob,
              ntot / njob, nchunks / njob))

    # run the stuff
    tmp = [s.splicables for s in segments]
    for s in segments: s.splicables = None  # poses not pickleable...
    with executor(max_workers=nworker) as pool:
        context = (sizes[end:], njob, segments, end, criteria, thresh)
        args = [range(njob)] + [it.repeat(context)]
        chunks = tqdm_parallel_map(
            pool, _grow_chunks, *args,
            unit='K worms', ascii=0, desc='growing worms',
            unit_scale=(ntot / njob / 1000))
        chunks = [x for x in chunks if x is not None]
    for s, t in zip(segments, tmp): s.splicables = t  # put the poses back

    # compose and sort results
    scores = np.concatenate([c[0] for c in chunks])
    order = np.argsort(scores)
    scores = scores[order]
    lowidx = np.concatenate([c[1] for c in chunks])[order]
    lowpos = np.concatenate([c[2] for c in chunks])[order]
    score_check = sum(c.score([lowpos[:, i]
                               for i in range(len(segments))]) for c in criteria)
    assert np.allclose(score_check, scores)
    return Worms(segments, scores, lowidx, lowpos)
