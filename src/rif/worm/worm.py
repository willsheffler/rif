from numpy.linalg import inv
from pyrosetta import rosetta as ros
from rif import rcl, homog, vis
from tqdm import tqdm
import functools as ft
import itertools as it
from concurrent.futures import ThreadPoolExecutor, as_completed
import multiprocessing
import os
import sys
import abc
import numpy as np
import multiprocessing

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


def trim_pose(pose, resid, direction, pad=0):
    "trim end of pose from direction, leaving <=pad residues beyond resid"
    if direction not in "NC":
        raise ValueError("direction must be 'N' or 'C'")
    if not 0 < resid <= len(pose):
        raise ValueError("resid %i out of bounds %i" % (resid, len(pose)))
    p = ros.core.pose.Pose()
    if direction == 'N':
        lb, ub = max(resid - pad, 1), len(pose)
    elif direction == 'C':
        lb, ub = 1, min(resid + pad, len(pose))
    # print('trim_pose lbub', lb, ub, 'len', len(pose), 'resid', resid)
    ros.core.pose.append_subpose_to_pose(p, pose, lb, ub)
    return p


def worst_CN_connect(p):
    for ir in range(1, len(p)):
        worst = 0
        if (p.residue(ir).is_protein() and
                p.residue(ir + 1).is_protein() and not (
                ros.core.pose.is_upper_terminus(p, ir) or
                ros.core.pose.is_lower_terminus(p, ir + 1))):
            dist = p.residue(ir).xyz('C').distance(p.residue(ir + 1).xyz('N'))
            worst = max(abs(dist - 1.32), worst)
    return worst


def reorder_spliced_as_N_to_C(body_chains, polarities):
    "remap chains of each body such that concatenated chains are N->C"
    if len(body_chains) != len(polarities) + 1:
        raise ValueError('must be one more body_chains than polarities')
    chains, pol = [[]], {}
    if not all(0 < len(dg) for dg in body_chains):
        raise ValueError('body_chains values must be [enterexit], '
                         '[enter,exit], or [enter, ..., exit')
    for i in range(1, len(polarities)):
        if len(body_chains[i]) == 1:
            if polarities[i - 1] != polarities[i]:
                raise ValueError('polarity mismatch on single chain connect')
    for i, dg in enumerate(body_chains):
        chains[-1].append(dg[0])
        if i != 0: pol[len(chains) - 1] = polarities[i - 1]
        if len(dg) > 1: chains.extend([x] for x in dg[1:])
    for i, chain in enumerate(chains):
        if i in pol and pol[i] == 'C':
            chains[i] = chains[i][::-1]
    return chains


class WormCriteria(abc.ABC):

    @abc.abstractmethod
    def score(self): pass

    def canonical_alignment(self): return

    def last_body_same_as(self): return

    def check_topolopy(self, segments): return


class AxesIntersect(WormCriteria):

    def __init__(self, *args, **kw):
        raise NotImplementedError


class SegmentSym(WormCriteria):

    def __init__(self, symmetry, from_seg=0, *, tol=1.0,
                 origin_seg=None, lever=100.0, to_seg=-1):
        self.symmetry = symmetry
        self.tol = tol
        self.from_seg = from_seg
        self.origin_seg = origin_seg
        self.lever = lever
        self.to_seg = to_seg
        self.rot_tol = tol / lever
        if self.symmetry[0] in 'cC':
            self.nfold = int(self.symmetry[1:])
            if self.nfold <= 0:
                raise ValueError('invalid symmetry: ' + symmetry)
            self.symangle = np.pi * 2.0 / self.nfold
        else: raise ValueError('can only do Cx symmetry for now')
        assert not origin_seg
        if self.tol <= 0: raise ValueError('tol should be > 0')

    def score(self, segpos, *, debug=False, **kwarge):
        x_from = segpos[self.from_seg]
        x_to = segpos[self.to_seg]
        xhat = x_to @ inv(x_from)
        trans = xhat[..., :, 3]
        if self.origin_seg:
            raise NotImplementedError
        elif self.nfold is 1:
            angle = homog.angle_of(xhat)
            carterrsq = np.sum(trans[..., :3]**2, axis=-1)
            roterrsq = angle**2
        else:
            axis, angle = homog.axis_angle_of(xhat)
            carterrsq = np.sum(trans * axis, axis=-1)**2
            roterrsq = (angle - self.symangle)**2
        return np.sqrt(carterrsq / self.tol**2 + roterrsq / self.rot_tol**2)

    def canonical_alignment(self, segpos, **kwargs):
        x_from = segpos[self.from_seg]
        x_to = segpos[self.to_seg]
        xhat = x_to @ inv(x_from)
        axis, ang, cen = homog.axis_ang_cen_of(xhat)
        dotz = homog.hdot(axis, [0, 0, 1])[..., None]
        tgtaxis = np.where(dotz > 0, [0, 0, 1, 0], [0, 0, -1, 0])
        align = homog.hrot((axis + tgtaxis) / 2, np.pi, cen)
        align[..., :3, 3] -= cen[..., :3]
        return align

    def last_body_same_as(self):
        return self.from_seg

    def check_topolopy(self, segments):
        "for cyclic, global entry can't be same as global exit"
        # todo: should check this...
        # fromseg = segments[self.from_seg]
        # toseg = segments[self.to_seg]
        return


class SpliceSite:

    def __init__(self, sele, polarity):
        if isinstance(sele, str) or isinstance(sele, int):
            sele = [sele]
        self.selections = list(sele)
        self.polarity = polarity

    def resid(self, id, pose):
        resid = id if id >= 0 else len(pose) + 1 + id
        if not 0 < resid <= len(pose):
            raise ValueError('resid ' + str(resid)
                             + ' invalid for pose of size '
                             + str(len(pose)))
        return resid

    def resids_impl(self, sele, spliceable):
        if isinstance(sele, int):
            return set([self.resid(sele, spliceable.body)])
        elif isinstance(sele, str):
            x = sele.split(',')
            s = x[-1].split(':')
            chain = int(x[0]) if len(x) == 2 else None
            pose = spliceable.chains[chain] if chain else spliceable.body
            start = self.resid(int(s[0] or 1), pose)
            stop = self.resid(int(s[1] or -1), pose)
            step = int(s[2]) if len(s) > 2 else 1
            # print(start, stop + 1, step)
            resids = set()
            for ir in range(start, stop + 1, step):
                assert 0 < ir <= len(pose)
                resids.add(spliceable.start_of_chain[chain] + ir)
            return resids
        elif sele is None:
            return set([None])
        else:
            raise ValueError('selection must be int, str, or None')

    def resids(self, spliceabe):
        resids = set()
        for sele in self.selections:
            try:
                resids |= self.resids_impl(sele, spliceabe)
            except ValueError as e:
                raise ValueError('Error with selection '
                                 + str(sele) + ': ' + str(e))
        resids = sorted(resids)
        if not resids:
            raise ValueError('empty SpliceSite')
        return resids

    def __repr__(ss):
        return 'SpliceSite(' + str(ss.selections) + ', ' + ss.polarity + ')'


class Spliceable:

    def __init__(self, body, sites, *, bodyid=None):
        self.body = body
        chains = list(body.split_by_chain())
        self.start_of_chain = {i + 1: sum(len(c) for c in chains[:i])
                               for i in range(len(chains))}
        self.start_of_chain[None] = 0
        self.chains = {i + 1: c for i, c in enumerate(chains)}
        self.bodyid = bodyid
        if callable(sites):
            sites = sites(body)
        if isinstance(sites, SpliceSite):
            sites = [sites]
        self.sites = list(sites)
        for i, site in enumerate(self.sites):
            if isinstance(site, str):
                raise ValueError('site currently must be (sele, polarity)')
            if not isinstance(site, SpliceSite):
                if not hasattr(site, '__iter__'):
                    self.sites[i] = (site,)
                self.sites[i] = SpliceSite(*site)

    def spliceable_positions(self):
        """selection of resids, and map 'global' index to selected index"""
        resid_subset = set()
        for site in self.sites:
            resid_subset |= set(site.resids(self))
        resid_subset = np.array(list(resid_subset))
        # really? must be an easier way to 'invert' a mapping in numpy?
        N = len(self.body) + 1
        val, idx = np.where(0 == (np.arange(N)[np.newaxis, :] -
                                  resid_subset[:, np.newaxis]))
        to_subset = np.array(N * [-1])
        to_subset[idx] = val
        assert np.all(to_subset[resid_subset] == np.arange(len(resid_subset)))
        return resid_subset, to_subset

    def __repr__(self):
        return ('Spliceable: body=(' + str(len(self.body)) + ',' +
                str(self.body).splitlines()[0].split('/')[-1] +
                '), sites=' + str([(s.resids(self), s.polarity) for s in self.sites]))


class Segment:

    def __init__(self, spliceables, entry=None, exit=None):
        self.entrypol = entry or None
        self.exitpol = exit or None
        if not spliceables:
            raise ValueError('spliceables must not be empty, spliceables =' +
                             str(spliceables))
        for s in spliceables:
            if not isinstance(s, Spliceable):
                raise ValueError('Segment can only accept list of Spliceable')
        self.init(spliceables, entry, exit)

    def __len__(self):
        return len(self.bodyid)

    def init(self, spliceables=None, entry=None, exit=None):
        if not (entry or exit):
            raise ValueError('at least one of entry/exit required')
        self.spliceables = list(spliceables) or self.spliceables
        self.entrypol = entry or self.entrypol or None
        self.exitpol = exit or self.exitpol or None
        # each array has all in/out pairs
        self.x2exit, self.x2orgn, self.bodyid = [], [], []
        self.entryresid, self.exitresid = [], []
        self.entrysiteid, self.exitsiteid = [], []
        # this whole loop is pretty inefficient, but that probably
        # doesn't matter much given the cost subsequent operations (?)
        for ibody, spliceable in enumerate(self.spliceables):
            resid_subset, to_subset = spliceable.spliceable_positions()
            bodyid = ibody if spliceable.bodyid is None else spliceable.bodyid
            # extract 'stubs' from body at selected positions
            # rif 'stubs' have 'extra' 'features'... the raw field is
            # just bog-standard homogeneous matrices
            stubs = rcl.bbstubs(spliceable.body, resid_subset)['raw']
            if len(resid_subset) != stubs.shape[0]:
                raise ValueError("no funny residues supported")
            stubs_inv = inv(stubs)
            entry_sites = (list(enumerate(spliceable.sites)) if self.entrypol else
                           [(-1, SpliceSite(sele=[None], polarity=self.entrypol))])
            exit_sites = (list(enumerate(spliceable.sites)) if self.exitpol else
                          [(-1, SpliceSite(sele=[None], polarity=self.exitpol))])
            for isite, entry_site in entry_sites:
                if entry_site.polarity == self.entrypol:
                    for jsite, exit_site in exit_sites:
                        if isite != jsite and exit_site.polarity == self.exitpol:
                            for ires in entry_site.resids(spliceable):
                                istub_inv = (identity44f4 if not ires
                                             else stubs_inv[to_subset[ires]])
                                ires = ires or -1
                                for jres in exit_site.resids(spliceable):
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

    def same_bodies_as(self, other):
        bodies1 = [s.body for s in self.spliceables]
        bodies2 = [s.body for s in other.spliceables]
        return bodies1 == bodies2

    def make_pose_chains(self, index, position=None, pad=(0, 0)):
        """returns (segchains, rest)
        segchains elems are [enterexitchain] or, [enterchain, ..., exitchain]
        rest holds other chains IFF entre and exit in same chain
        """
        spliceable = self.spliceables[self.bodyid[index]]
        pose = spliceable.body
        chains = spliceable.chains
        ir_en, ir_ex = self.entryresid[index], self.exitresid[index]
        ch_en = pose.chain(ir_en) if ir_en > 0 else None
        ch_ex = pose.chain(ir_ex) if ir_ex > 0 else None
        pl_en, pl_ex = self.entrypol, self.exitpol
        if ch_en: ir_en -= spliceable.start_of_chain[ch_en]
        if ch_ex: ir_ex -= spliceable.start_of_chain[ch_ex]
        assert ch_en or ch_ex
        rest = list(chains.values())
        if ch_en: rest.remove(chains[ch_en])
        if ch_en == ch_ex:
            assert len(rest) + 1 == len(chains)
            p = trim_pose(chains[ch_en], ir_en, self.entrypol, pad[0])
            iexit1 = ir_ex - (pl_ex == 'C') * (len(chains[ch_en]) - len(p))
            p = trim_pose(p, iexit1, pl_ex, pad[1] - 1)
            enex, rest = [p], rest
        else:
            if ch_ex: rest.remove(chains[ch_ex])
            p_en = [chains[ch_en]] if ch_en else []
            p_ex = [chains[ch_ex]] if ch_ex else []
            if p_en: p_en[0] = trim_pose(p_en[0], ir_en, self.entrypol, pad[0])
            if p_ex: p_ex[0] = trim_pose(p_ex[0], ir_ex, pl_ex, pad[1] - 1)
            enex, rest = p_en + rest + p_ex, []
        if position is not None:
            position = rcl.to_rosetta_stub(position)
            enex, rest = [p.clone() for p in enex], [p.clone() for p in rest]
            for p in it.chain(enex, rest):
                ros.protocols.sic_dock.xform_pose(p, position)
        return enex, rest


class Worms:

    def __init__(self, segments, scores, indices, positions, criteria):
        self.segments = segments
        self.scores = scores
        self.indices = indices
        self.positions = positions
        self.criteria = criteria

    def __init___(self):
        return len(self.scores)

    def pose(self, which, align=True, withend=True, join=True):
        "makes a pose for the ith worm"
        if hasattr(which, '__iter__'):
            return (self.pose(w) for w in which)
        print("Will needs to fix bb O/H position!")
        rm_lower_t = ros.core.pose.remove_lower_terminus_type_from_pose_residue
        rm_upper_t = ros.core.pose.remove_upper_terminus_type_from_pose_residue
        entryexits = [seg.make_pose_chains(self.indices[which][iseg],
                                           self.positions[which][iseg])
                      for iseg, seg in enumerate(self.segments)]
        entryexits, rest = zip(*entryexits)
        chainslist = reorder_spliced_as_N_to_C(
            entryexits, [s.entrypol for s in self.segments[1:]])
        pose = ros.core.pose.Pose()
        for chains in chainslist:
            ros.core.pose.append_pose_to_pose(pose, chains[0], True)
            for chain in chains[1:]:
                rm_upper_t(pose, len(pose))
                rm_lower_t(chain, 1)
                ros.core.pose.append_pose_to_pose(pose, chain, not join)
        for chain in it.chain(*rest):
            ros.core.pose.append_pose_to_pose(pose, chain, True)
        if align:
            align = [c.canonical_alignment(self.positions[which])
                     for c in self.criteria]
            align = [x for x in align if x is not None]
            assert len(align) < 2  # should this be allowed?
            if len(align) == 1:
                x = rcl.to_rosetta_stub(align[0])
                ros.protocols.sic_dock.xform_pose(pose, x)
        assert worst_CN_connect(pose) < 0.5
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


def _grow_chunk(samp, segpos, conpos, segs, end, criteria, thresh, matchlast):
    if matchlast is not None:
        ndimchunk = segpos[0].ndim - 2
        if matchlast < ndimchunk:
            bidA = segs[matchlast].bodyid
            bidB = segs[-1].bodyid[samp[-1]]
            idx = (slice(None),) * matchlast + (bidA == bidB,)
            segpos = segpos[: matchlast] + [x[idx] for x in segpos[matchlast:]]
            conpos = conpos[: matchlast] + [x[idx] for x in conpos[matchlast:]]
            idxmap = np.where(bidA == bidB)[0]
        elif segs[matchlast].bodyid[samp[matchlast - ndimchunk]] != bidB:
            return  # last body doesn't match for this whole chunk
    segpos, conpos = segpos[:end], conpos[:end]
    for iseg, seg in enumerate(segs[end:]):
        segpos.append(conpos[-1] @ seg.x2orgn[samp[iseg]])
        if seg is not segs[-1]:
            conpos.append(conpos[-1] @ seg.x2exit[samp[iseg]])
    score = sum(c.score(segpos=segpos) for c in criteria)
    ilow0 = np.where(score < thresh)
    sampidx = tuple(np.repeat(i, len(ilow0[0])) for i in samp)
    lowpostmp = []
    for iseg in range(len(segpos)):
        ilow = ilow0[: iseg + 1] + (0,) * (segpos[0].ndim - 2 - (iseg + 1))
        lowpostmp.append(segpos[iseg][ilow])
    ilow1 = (ilow0 if matchlast is None else ilow0[: matchlast] +
             (idxmap[ilow0[matchlast]],) + ilow0[matchlast + 1:])
    return score[ilow0], np.array(ilow1 + sampidx).T, np.stack(lowpostmp, 1)


def _grow_chunks(ijob, context):
    sampsizes, njob, segments, end, criteria, thresh, matchlast = context
    samples = it.product(*(range(n) for n in sampsizes))
    segpos, connpos = _chain_xforms(segments[:end])  # common data
    args = [list(samples)[ijob::njob]] + [it.repeat(x) for x in (
        segpos, connpos, segments, end, criteria, thresh, matchlast)]
    chunk = list(map(_grow_chunk, *args))
    return [np.concatenate([c[i] for c in chunk])
            for i in range(3)] if chunk else None


def grow(segments, criteria, *, thresh=2, expert=0, memlim=1e6,
         executor=None, max_workers=None, debug=0, jobmult=32, verbose=0):
    if verbose:
        print('grow')
        for i, seg in enumerate(segments):
            print(' segment', i, 'enter:', seg.entrypol, 'exit:', seg.exitpol)
            for sp in seg.spliceables:
                print('   ', sp)
    criteria = [criteria] if isinstance(criteria, WormCriteria) else criteria
    # checks and setup
    if segments[0].entrypol is not None:
        raise ValueError('beginning of worm cant have entry')
    if segments[-1].exitpol is not None:
        raise ValueError('end of worm cant have exit')
    for a, b in zip(segments[:-1], segments[1:]):
        if not (a.exitpol and b.entrypol and a.exitpol != b.entrypol):
            raise ValueError('incompatible exit->entry polarity: '
                             + str(a.exitpol) + '->'
                             + str(b.entrypol) + ' on segment pair: '
                             + str((segments.index(a), segments.index(b))))
    matchlast = ([c.last_body_same_as() for c in criteria] or [None])[0]
    if matchlast is not None and not expert and (
            not segments[matchlast].same_bodies_as(segments[-1])):
        raise ValueError("segments[matchlast] not same as segments[-1], "
                         + "if you're sure, pass expert=True")
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
          'worm/job = {:,}, chunk/job = {}, sizes={}'.format(
              ntot, nchunk, nchunks, nworker, njob,
              ntot / njob, nchunks / njob, sizes))

    # run the stuff
    tmp = [s.spliceables for s in segments]
    for s in segments: s.spliceables = None  # poses not pickleable...
    with executor(max_workers=nworker) as pool:
        context = (sizes[end:], njob, segments, end, criteria, thresh,
                   matchlast)
        args = [range(njob)] + [it.repeat(context)]
        chunks = tqdm_parallel_map(
            pool, _grow_chunks, *args,
            unit='K worms', ascii=0, desc='growing worms',
            unit_scale=int(ntot / njob / 1000))
        chunks = [x for x in chunks if x is not None]
    for s, t in zip(segments, tmp): s.spliceables = t  # put the poses back

    # compose and sort results
    scores = np.concatenate([c[0] for c in chunks])
    order = np.argsort(scores)
    scores = scores[order]
    lowidx = np.concatenate([c[1] for c in chunks])[order]
    lowpos = np.concatenate([c[2] for c in chunks])[order]
    score_check = sum(c.score([lowpos[:, i]
                               for i in range(len(segments))]) for c in criteria)
    assert np.allclose(score_check, scores)
    return Worms(segments, scores, lowidx, lowpos, criteria)
