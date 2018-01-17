from numpy.linalg import inv
from pyrosetta import rosetta as ros
from pyrosetta.rosetta.core import scoring
from rif import rcl, homog, vis, sym
from tqdm import tqdm  # progress bar utility
import functools as ft
import itertools as it
from concurrent.futures import ThreadPoolExecutor, as_completed
import multiprocessing
import os
import sys
import abc
import numpy as np
import multiprocessing


identity44f8 = np.identity(4, dtype='f4')
identity44f8 = np.identity(4, dtype='f8')
Ux = np.array([1, 0, 0, 0])
Uy = np.array([0, 1, 0, 0])
Uz = np.array([0, 0, 1, 0])

# todo: the following should go elsewhere...


def cpu_count():
    try: return int(os.environ['SLURM_CPUS_ON_NODE'])
    except: return multiprocessing.cpu_count()


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


def symfile_path(name):
    path, _ = os.path.split(__file__)
    return os.path.join(path, 'rosetta_symdef', name + '.sym')


class WormCriteria(abc.ABC):

    @abc.abstractmethod
    def score(self): pass

    def canonical_alignment(self, *args, **kw): return identity44f8

    def last_body_same_as(self): return

    def check_topolopy(self, segments): return

    def is_cyclic(self): return False

    def sym_axes(self): return []

    def sym_frames(self): return [identity44f8]

    def symdef(self): return None


class AxesIntersect(WormCriteria):

    def __init__(self, symname, tgtaxis1, tgtaxis2, from_seg, *, tol=1.0,
                 lever=50, to_seg=-1, frames=None, distinct_axes=False):
        self.symname = symname
        self.from_seg = from_seg
        if len(tgtaxis1) == 2: tgtaxis1 += [0, 0, 0, 1],
        if len(tgtaxis2) == 2: tgtaxis2 += [0, 0, 0, 1],
        self.tgtaxis1 = (tgtaxis1[0], homog.hnormalized(tgtaxis1[1]),
                         homog.hpoint(tgtaxis1[2]))
        self.tgtaxis2 = (tgtaxis2[0], homog.hnormalized(tgtaxis2[1]),
                         homog.hpoint(tgtaxis2[2]))
        assert 3 == len(self.tgtaxis1)
        assert 3 == len(self.tgtaxis2)
        self.angle = homog.angle(tgtaxis1[1], tgtaxis2[1])
        self.tol = tol
        self.lever = lever
        self.to_seg = to_seg
        self.rot_tol = tol / lever
        self.frames = frames
        if frames is None:
            x1 = homog.hrot(tgtaxis1[1], 360 / tgtaxis1[0], self.tgtaxis1[2])
            x2 = homog.hrot(tgtaxis2[1], 360 / tgtaxis2[0], self.tgtaxis2[2])
            self.frames = [identity44f8, x1, x2, x1 @ x2]
        self.distinct_axes = distinct_axes  # -z not same as z (for T33)

    def score(self, segpos, verbose=False, **kw):
        cen1 = segpos[self.from_seg][..., :, 3]
        cen2 = segpos[self.to_seg][..., :, 3]
        ax1 = segpos[self.from_seg][..., :, 2]
        ax2 = segpos[self.to_seg][..., :, 2]
        if self.distinct_axes:
            p, q = homog.line_line_closest_points_pa(cen1, ax1, cen2, ax2)
            dist = homog.hnorm(p - q)
            cen = (p + q) / 2
            ax1c = homog.hnormalized(cen1 - cen)
            ax2c = homog.hnormalized(cen2 - cen)
            ax1 = np.where(homog.hdot(ax1, ax1c)[..., None] > 0, ax1, -ax1)
            ax2 = np.where(homog.hdot(ax2, ax2c)[..., None] > 0, ax2, -ax2)
            ang = np.arccos(homog.hdot(ax1, ax2))
        else:
            dist = homog.line_line_distance_pa(cen1, ax1, cen2, ax2)
            ang = np.arccos(np.abs(homog.hdot(ax1, ax2)))
        roterr2 = (ang - self.angle)**2
        return np.sqrt(roterr2 / self.rot_tol**2 + (dist / self.tol)**2)

    def canonical_alignment(self, segpos, debug=0, **kw):
        cen1 = segpos[self.from_seg][..., :, 3]
        cen2 = segpos[self.to_seg][..., :, 3]
        ax1 = segpos[self.from_seg][..., :, 2]
        ax2 = segpos[self.to_seg][..., :, 2]
        if not self.distinct_axes and homog.angle(ax1, ax2) > np.pi / 2:
            ax2 = -ax2
        p, q = homog.line_line_closest_points_pa(cen1, ax1, cen2, ax2)
        cen = (p + q) / 2
        # ax1 = homog.hnormalized(cen1 - cen)
        # ax2 = homog.hnormalized(cen2 - cen)
        x = homog.align_vectors(ax1, ax2, self.tgtaxis1[1], self.tgtaxis2[1])
        x[..., :, 3] = - x @cen
        if debug:
            print('angs', homog.angle_degrees(ax1, ax2),
                  homog.angle_degrees(self.tgtaxis1[1], self.tgtaxis2[1]))
            print('ax1', ax1)
            print('ax2', ax2)
            print('xax1', x @ ax1)
            print('tax1', self.tgtaxis1[1])
            print('xax2', x @ ax2)
            print('tax2', self.tgtaxis2[1])
            raise AssertionError
            # if not (np.allclose(x @ ax1, self.tgtaxis1[1], atol=1e-2) and
            #         np.allclose(x @ ax2, self.tgtaxis2[1], atol=1e-2)):
            #     print(homog.angle(self.tgtaxis1[1], self.tgtaxis2[1]))
            #     print(homog.angle(ax1, ax2))
            #     print(x @ ax1)
            #     print(self.tgtaxis1[1])
            #     print(x @ ax2)
            #     print(self.tgtaxis2[1])
            #     raise AssertionError('homog.align_vectors sucks')

        return x

    def sym_frames(self): return self.frames

    def symdef(self): return symfile_path(self.symname)

    def sym_axes(self): return [self.tgtaxis1, self.tgtaxis2]


def D2(c2=0, c2b=-1, **kw):
    return AxesIntersect('D2', (2, Uz), (2, Ux), c2, to_seg=c2b, **kw)


def D3(c3=0, c2=-1, **kw):
    return AxesIntersect('D3', (3, Uz), (2, Ux), c3, to_seg=c2, **kw)


def D4(c4=0, c2=-1, **kw):
    return AxesIntersect('D4', (4, Uz), (2, Ux), c4, to_seg=c2, **kw)


def D5(c5=0, c2=-1, **kw):
    return AxesIntersect('D5', (5, Uz), (2, Ux), c5, to_seg=c2, **kw)


def D6(c6=0, c2=-1, **kw):
    return AxesIntersect('D6', (6, Uz), (2, Ux), c6, to_seg=c2, **kw)


def Tetrahedral(c3=None, c2=None, c3b=None, **kw):
    if 1 is not (c3b is None) + (c3 is None) + (c2 is None):
        raise ValueError('must specify exactly two of c3, c2, c3b')
    if c2 is None: from_seg, to_seg, nf1, nf2, ex = c3b, c3, 7, 3, 2
    if c3 is None: from_seg, to_seg, nf1, nf2, ex = c3b, c2, 7, 2, 3
    if c3b is None: from_seg, to_seg, nf1, nf2, ex = c3, c2, 3, 2, 7
    return AxesIntersect('T', from_seg=from_seg, to_seg=to_seg,
                         tgtaxis1=(max(3, nf1), sym.tetrahedral_axes[nf1]),
                         tgtaxis2=(max(3, nf2), sym.tetrahedral_axes[nf2]),
                         frames=sym.tetrahedral_frames,
                         distinct_axes=(nf1 == 7), **kw)


def Octahedral(c4=None, c3=None, c2=None, **kw):
    if 1 is not (c4 is None) + (c3 is None) + (c2 is None):
        raise ValueError('must specify exactly two of c4, c3, c2')
    if c2 is None: from_seg, to_seg, nf1, nf2, ex = c4, c3, 4, 3, 2
    if c3 is None: from_seg, to_seg, nf1, nf2, ex = c4, c2, 4, 2, 3
    if c4 is None: from_seg, to_seg, nf1, nf2, ex = c3, c2, 3, 2, 4
    return AxesIntersect('O', from_seg=from_seg, to_seg=to_seg,
                         tgtaxis1=(nf1, sym.octahedral_axes[nf1]),
                         tgtaxis2=(nf2, sym.octahedral_axes[nf2]),
                         frames=sym.octahedral_frames, **kw)


def Icosahedral(c5=None, c3=None, c2=None, **kw):
    if 1 is not (c5 is None) + (c3 is None) + (c2 is None):
        raise ValueError('must specify exactly two of c5, c3, c2')
    if c2 is None: from_seg, to_seg, nf1, nf2, ex = c5, c3, 5, 3, 2
    if c3 is None: from_seg, to_seg, nf1, nf2, ex = c5, c2, 4, 2, 3
    if c5 is None: from_seg, to_seg, nf1, nf2, ex = c3, c2, 3, 2, 5
    return AxesIntersect('I', from_seg=from_seg, to_seg=to_seg,
                         tgtaxis1=(nf1, sym.icosahedral_axes[nf1]),
                         tgtaxis2=(nf2, sym.icosahedral_axes[nf2]),
                         frames=sym.icosahedral_frames, **kw)


class Cyclic(WormCriteria):

    def __init__(self, symmetry, from_seg=0, *, tol=1.0,
                 origin_seg=None, lever=50.0, to_seg=-1,
                 relweight=1.0):
        self.symmetry = symmetry
        self.tol = tol
        self.from_seg = from_seg
        self.origin_seg = origin_seg
        self.lever = lever
        self.to_seg = to_seg
        self.rot_tol = tol / lever
        self.relweight = relweight if abs(relweight) > 0.001 else None
        if self.symmetry[0] in 'cC':
            self.nfold = int(self.symmetry[1:])
            if self.nfold <= 0:
                raise ValueError('invalid symmetry: ' + symmetry)
            self.symangle = np.pi * 2.0 / self.nfold
        else: raise ValueError('can only do Cx symmetry for now')
        assert not origin_seg
        if self.tol <= 0: raise ValueError('tol should be > 0')
        self.frames = [homog.hrot(Uz, self.symangle * i)
                       for i in range(self.nfold)]

    def score(self, segpos, *, verbose=False, **kw):
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
            carterrsq = homog.hdot(trans, axis)**2
            if self.relweight is not None:
                distsq = np.sum(trans[..., :3]**2, axis=-1)
                relerrsq = carterrsq / distsq
                relerrsq[np.isnan(relerrsq)] = 9e9
                carterrsq += self.relweight * relerrsq  # too much of a hack??
            if verbose:
                print('axis', axis[0])
                print('trans', trans[0])
                print('dot trans', homog.hdot(trans, axis)[0])
                print('ang', angle[0] * 180 / np.pi)
            roterrsq = (angle - self.symangle)**2
            # if verbose:
            # print('cart', carterrsq)
            # print('rote', roterrsq)
            # print('ang', angle * 180 / np.pi)
        return np.sqrt(carterrsq / self.tol**2 +
                       roterrsq / self.rot_tol**2)

    def canonical_alignment(self, segpos, **kwargs):
        x_from = segpos[self.from_seg]
        x_to = segpos[self.to_seg]
        xhat = x_to @ inv(x_from)
        axis, ang, cen = homog.axis_ang_cen_of(xhat)
        # print('aln', axis)
        # print('aln', ang * 180 / np.pi)
        # print('aln', cen)
        # print('aln', xhat[..., :, 3])
        dotz = homog.hdot(axis, Uz)[..., None]
        tgtaxis = np.where(dotz > 0, [0, 0, 1, 0], [0, 0, -1, 0])
        align = homog.hrot((axis + tgtaxis) / 2, np.pi, cen)
        align[..., :3, 3] -= cen[..., :3]
        return align

    def last_body_same_as(self): return self.from_seg

    def check_topolopy(self, segments):
        "for cyclic, global entry can't be same as global exit"
        # todo: should check this...
        # fromseg = segments[self.from_seg]
        # toseg = segments[self.to_seg]
        return

    def is_cyclic(self): return True

    def sym_frames(self): return self.frames

    def symdef(self): return symfile_path('C' + str(self.nfold))

    def sym_axes(self): return [(self.nfold, Uz, [0, 0, 0, 1])]


class SpliceSite:

    def __init__(self, sele, polarity, chain=None):
        if isinstance(sele, str) or isinstance(sele, int):
            sele = [sele]
        self.selections = list(sele)
        self.polarity = polarity
        self.chain = chain

    def resid(self, id, pose):
        resid = id if id >= 0 else len(pose) + 1 + id
        if not 0 < resid <= len(pose):
            raise ValueError('resid ' + str(resid)
                             + ' invalid for pose of size '
                             + str(len(pose)))
        return resid

    def resids_impl(self, sele, spliceable):
        if isinstance(sele, int):
            if self.chain is None:
                return set([self.resid(sele, spliceable.body)])
            else:
                ir = self.resid(sele, spliceable.chains[self.chain])
                ir += spliceable.start_of_chain[self.chain]
                return set([ir])
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

    def __repr__(self):
        c = '' if self.chain is None else ', chain=' + str(self.chain)
        return 'SpliceSite(' + str(self.selections) + \
            ', ' + self.polarity + c + ')'


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
                if isinstance(site, dict):
                    self.sites[i] = SpliceSite(**site)
                else:
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
            raise ValueError('spliceables must not be empty, spliceables ='
                             + str(spliceables))
        for s in spliceables:
            if not isinstance(s, Spliceable):
                raise ValueError('Segment only accepts list of Spliceable')
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
            stubs = stubs.astype('f8')
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
                                istub_inv = (identity44f8 if not ires
                                             else stubs_inv[to_subset[ires]])
                                ires = ires or -1
                                for jres in exit_site.resids(spliceable):
                                    jstub = (identity44f8 if not jres
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

    def make_pose_chains(self, index, position=None, pad=(0, 0), ):
        """returns (segchains, rest)
        segchains elems are [enterexitchain] or, [enterchain, ..., exitchain]
        rest holds other chains IFF enter and exit in same chain
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


def _subpose(pose, lb, ub):
    assert lb > 0 and ub <= len(pose)
    p = ros.core.pose.Pose()
    ros.core.pose.append_subpose_to_pose(p, pose, lb, ub)
    return p


def _cyclic_permute_chains(chainslist, polarity, spliceres):
    rm_lower_t = ros.core.pose.remove_lower_terminus_type_from_pose_residue
    rm_upper_t = ros.core.pose.remove_upper_terminus_type_from_pose_residue
    beg, end = chainslist[0], chainslist[-1]
    n2c = polarity == 'N'
    if n2c:
        stub1 = rcl.bbstubs(beg[0], [spliceres])
        stub2 = rcl.bbstubs(end[-1], [len(end[-1])])
        beg[0] = _subpose(beg[0], spliceres + 1, len(beg[0]))
        rm_lower_t(beg[0], 1)
        rm_upper_t(end[-1], len(end[-1]))
    else:
        stub1 = rcl.bbstubs(beg[-1], [spliceres])
        stub2 = rcl.bbstubs(end[0], [1])
        beg[-1] = _subpose(beg[-1], 1, spliceres - 1)
        rm_lower_t(beg[-1], len(beg[-1]))
        rm_upper_t(end[0], 1)
    xalign = stub1['raw'][0] @ np.linalg.inv(stub2['raw'][0])
    for p in end: rcl.xform_pose(xalign, p)
    if n2c: chainslist[0] = end + beg
    else: chainslist[0] = beg + end
    chainslist = chainslist[:-1]


class Worms:

    def __init__(self, segments, scores, indices, positions, criteria, detail):
        self.segments = segments
        self.scores = scores
        self.indices = indices
        self.positions = positions
        self.criteria = criteria
        self.detail = detail
        self.score0 = scoring.symmetry.symmetrize_scorefunction(
            scoring.ScoreFunctionFactory.create_score_function('score0'))
        self.splicepoint_cache = {}

    def pose(self, which, *, align=True, end=None, join=True,
             only_connected=None, cyclic_permute=False, **kw):
        "makes a pose for the ith worm"
        if hasattr(which, '__iter__'): return (
            self.pose(w, align=align, end=end, join=join,
                      only_connected=only_connected, **kw) for w in which)
        # print("Will needs to fix bb O/H position!")
        rm_lower_t = ros.core.pose.remove_lower_terminus_type_from_pose_residue
        rm_upper_t = ros.core.pose.remove_upper_terminus_type_from_pose_residue
        if end is None:
            end = not self.criteria[0].is_cyclic()
        if only_connected is None:
            only_connected = not self.criteria[0].is_cyclic()
        if cyclic_permute is None:
            cyclic_permute = self.criteria[0].is_cyclic()
        elif cyclic_permute and not self.criteria[0].is_cyclic():
            raise ValueError('cyclic_permute should only be used for Cyclic')
        iend = None if end else -1
        entryexits = [seg.make_pose_chains(self.indices[which][iseg],
                                           self.positions[which][iseg])
                      for iseg, seg in enumerate(self.segments[:iend])]
        entryexits, rest = zip(*entryexits)
        chainslist = reorder_spliced_as_N_to_C(
            entryexits, [s.entrypol for s in self.segments[1:iend]])
        if align:
            align = [c.canonical_alignment(self.positions[which], **kw)
                     for c in self.criteria]
            align = [x for x in align if x is not None]
            assert len(align) < 2  # should this be allowed?
            if len(align) == 1:
                for p in it.chain(*chainslist, *rest):
                    rcl.xform_pose(align[0], p)
        if cyclic_permute and len(chainslist) > 1:
            # todo: this is only correct if 1st seg is one chain
            spliceres = self.segments[-1].entryresid[self.indices[which, -1]]
            bodyid = self.segments[0].bodyid[self.indices[which, 0]]
            origlen = len(self.segments[0].spliceables[bodyid].body)
            i = -1 if self.segments[-1].entrypol == 'C' else 0
            spliceres -= origlen - len(chainslist[0][i])
            _cyclic_permute_chains(chainslist, self.segments[-1].entrypol,
                                   spliceres + 1)
        pose = ros.core.pose.Pose()
        splicepoints = []
        for chains in chainslist:
            if only_connected and len(chains) is 1: continue
            ros.core.pose.append_pose_to_pose(pose, chains[0], True)
            for chain in chains[1:]:
                rm_upper_t(pose, len(pose))
                rm_lower_t(chain, 1)
                splicepoints.append(len(pose))
                ros.core.pose.append_pose_to_pose(pose, chain, not join)
        self.splicepoint_cache[which] = splicepoints
        if not only_connected:
            for chain in it.chain(*rest):
                ros.core.pose.append_pose_to_pose(pose, chain, True)
        assert rcl.worst_CN_connect(pose) < 0.5
        return pose

    def splicepoints(self, which):
        if not which in self.splicepoint_cache:
            self.pose(which)
        return self.splicepoint_cache[which]

    def sympose(self, which, fullatom=False, score=False):
        if hasattr(which, '__iter__'): return (self.sympose(w) for w in which)
        p = self.pose(which, splicepoints=True)
        if not fullatom:
            ros.core.util.switch_to_residue_type_set(p, 'centroid')
        ros.core.pose.symmetry.make_symmetric_pose(
            p, self.criteria[0].symdef())
        if not score: return p
        return p, self.score0(p)

    def __len__(self): return len(self.scores)

    def __getitem__(self, i):
        return (i, self.scores[i],) + self.sympose(i, score=True)


def _chain_xforms(segments):
    os.environ['OMP_NUM_THREADS'] = '1'
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
    os.environ['OMP_NUM_THREADS'] = '1'
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
    os.environ['OMP_NUM_THREADS'] = '1'
    sampsizes, njob, segments, end, criteria, thresh, matchlast = context
    samples = it.product(*(range(n) for n in sampsizes))
    segpos, connpos = _chain_xforms(segments[:end])  # common data
    args = [list(samples)[ijob::njob]] + [it.repeat(x) for x in (
        segpos, connpos, segments, end, criteria, thresh, matchlast)]
    chunk = list(map(_grow_chunk, *args))
    return [np.concatenate([c[i] for c in chunk])
            for i in range(3)] if chunk else None


def grow(segments, criteria, *, thresh=2, expert=0, memsize=1e6,
         executor=None, max_workers=None, verbose=0, jobmult=32,
         chunklim=None):

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
        executor = ThreadPoolExecutor  # todo: some kind of null executor?
        max_workers = 1
    if max_workers is None: max_workers = cpu_count()
    sizes = [len(s.bodyid) for s in segments]
    end = len(segments) - 1
    while end > 1 and (np.prod(sizes[end:]) < max_workers or
                       memsize <= 64 * np.prod(sizes[:end])): end -= 1
    ntot, chunksize, nchunks = (np.product(x)
                                for x in (sizes, sizes[:end], sizes[end:]))
    nworker = max_workers or cpu_count()
    njob = nworker * jobmult
    njob = min(njob, nchunks)
    print('tot = {:,}, chunksize = {:,}, nchunks = {:,}, nworker = {}, '
          'njob = {}, worm/job = {:,}, chunk/job = {}, sizes={}'.format(
              ntot, chunksize, nchunks, nworker, njob,
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
    lowposlist = [lowpos[:, i] for i in range(len(segments))]
    score_check = sum(c.score(lowposlist, verbose=verbose) for c in criteria)
    assert np.allclose(score_check, scores)
    detail = dict(ntot=ntot, chunksize=chunksize, nchunks=nchunks,
                  nworker=nworker, njob=njob, sizes=sizes, end=end)
    return Worms(segments, scores, lowidx, lowpos, criteria, detail)
