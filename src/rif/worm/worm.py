


class SplicePoint:
    def __init__(self, residue_numbers, direction):
        self.residue_numbers = residue_numbers
        self.direction = direction


class Segment:
    def __init__(self, pose, splicepoints):
        self.pose = pose
        self.splices = splicepoints


class SegmentSet:
    def __init__(self, segments):
        self.segments = segments


class SegmentSets:
    def __init__(segment_sets):
        self.segment_sets = segment_sets


class XformConstraint:
    def __init__(self, from_segment=0, to_segment=-1, sym='C1',
                has_same_sym_as=None, is_exactly=None, tol=0.1, lever=100.0):
        self.from_segment = from_segment
        self.to_segment = to_segment
        self.sym = sym
        self.has_same_sym_as = has_same_sym_as
        self.is_exactly = is_exactly
        self.tol = tol
        self.lever = lever


def grow_worms(segment_sets, constraints, score_funcs=None):
    return None