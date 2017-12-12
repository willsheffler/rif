from rif import rcl
from rif.eigen_types import x3_identity, x3_inverse
import numpy as np


class SplicePositions:

    def __init__(self, residue_numbers, direction):
        self.residue_numbers = list(residue_numbers)
        self.direction = direction

    def __len__(self):
        return len(self.residue_numbers)

    def __iter__(self):
        return iter(self.residue_numbers)


class Segment:

    def __init__(self, pose, splice_groups, entry_dir=None, exit_dir=None):
        assert entry_dir or exit_dir
        self.pose = pose
        self.splice_groups = list(splice_groups)
        self.entry_dir = entry_dir
        self.exit_dir = exit_dir
        self.update_geometry()

    def splice_position_map(self):
        """selection of resis, and map 'global' index to selected index"""
        resi_sele = set()
        for splice_positions in self.splice_groups:
            resi_sele |= set(splice_positions)
        resi_sele = np.array(list(resi_sele))
        # really? must be an easier way to 'invert' a selection in numpy
        idx_range = np.arange(self.pose.size() + 1)
        val, idx = np.where((np.expand_dims(idx_range, 0) -
                             np.expand_dims(resi_sele, 1)) == 0)
        sele_inv = np.array((self.pose.size() + 1) * [-1])
        sele_inv[idx] = val
        assert (sele_inv[resi_sele] == np.arange(len(resi_sele))).all()
        return resi_sele, sele_inv

    def update_geometry(self):
        resi_sele, sele_inv = self.splice_position_map()
        # extract 'stubs' from pose at selected positions
        # rif 'stubs' have 'extra' 'features'... the raw field is
        # just bog-standard homogeneous matrices
        bbstubs = rcl.bbstubs(self.pose, resi_sele)['raw']
        assert len(resi_sele) is bbstubs.shape[
            0], "no funny residues supported"
        bbstubs_inv = np.linalg.inv(bbstubs)
        self.entry2exit = list()  # all in/out xforms
        self.entry2cntr = list()  # all in to centry (redundant)
        self.entry_resi = list()
        self.exit_resi = list()
        if self.entry_dir is None:  # beginning of 'worm'
            for group in self.splice_groups:
                if group.direction == self.exit_dir:
                    for jr in group:
                        self.entry2exit.append(bbstubs[sele_inv[jr]])
                        self.entry2cntr.append(x3_identity)
                        self.entry_resi.append(np.nan)
                        self.exit_resi.append(jr)
        elif self.exit_dir is None:  # at end of 'worm'
            for group in self.splice_groups:
                if group.direction == self.entry_dir:
                    for ir in group:
                        self.entry2exit.append(bbstubs_inv[sele_inv[ir]])
                        self.entry2cntr.append(bbstubs_inv[sele_inv[ir]])
                        self.entry_resi.append(ir)
                        self.exit_resi.append(np.nan)
        else:  # in middle of worm
            for i, group1 in enumerate(self.splice_groups):
                if group1.direction != self.entry_dir:
                    continue
                for j, group2 in enumerate(self.splice_groups):
                    if group2.direction != self.exit_dir:
                        continue
                    if i == j:
                        continue
                    for ir in group1:
                        print(ir)
                        print(sele_inv[ir])
                        istub = bbstubs_inv[sele_inv[ir]]
                        for jr in group2:
                            jstub = bbstubs[sele_inv[jr]]
                            self.entry2exit.append(istub @ jstub)
                            self.entry2cntr.append(istub)
                            self.entry_resi.append(ir)
                            self.exit_resi.append(jr)
        assert len(self.entry2exit) > 0
        self.entry2exit = np.stack(self.entry2exit)
        self.entry2cntr = np.stack(self.entry2cntr)
        self.entry_resi = np.array(self.entry_resi)
        self.exit_resi = np.array(self.exit_resi)


class SegmentSet:

    def __init__(self, segments):
        assert len(segments) > 0
        self.entry_dir = segments[0].entry_dir
        self.exit_dir = segments[0].exit_dir
        for s in segments:
            assert s.entry_dir == self.entry_dir
            assert s.exit_dir == self.exit_dir
        self.segments = segments


class SegmentSets:

    def __init__(segment_sets):
        self.segment_sets = segment_sets


class TransformConstraint:

    def __init__(self, from_segment=0, to_segment=-1, sym='C1',
                 has_same_sym_as=None, is_exactly=None, tol=0.1, lever=100.0):
        self.from_segment = from_segment
        self.to_segment = to_segment
        self.sym = sym
        self.has_same_sym_as = has_same_sym_as
        self.is_exactly = is_exactly
        self.tol = tol
        self.lever = lever

    def is_satisfactory(positions):
        return False


def grow_worms(segment_sets, constraints, score_funcs=None, cache=None):
    return None
