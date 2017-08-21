from rif import stripe_index_3d, rcl
if rcl.HAVE_PYROSETTA:
    from rosetta.numeric import xyzVector_float_t
    from rosetta.core.kinematics import Stub


class PoseResnumIndex(object):
    def __init__(self, pose, clash_dis=3.2, nbr_dis=7.0, nbr_atom='CA'):
        bb_atoms = rcl.atoms(pose, 'bbheavy')
        self.clash_check = stripe_index_3d(clash_dis, bb_atoms)
        cens, resnums = list(), list()
        for res in pose:
            if res.is_protein():
                cens.append(res.xyz('CA'))
                resnums.append(res.seqpos())
        self.resi_map = stripe_index_3d(nbr_dis, cens, resnums)
        self.nbr_atom = nbr_atom

    def clashes(self, pose):
        for res in pose:
            for ia in range(1, res.nheavyatoms() + 1):
                xyz = res.xyz(ia)
                if self.clash_check.neighbor_exists(xyz):
                    return True
        return False

    def contacting_resnums(self, pose):
        if self.clashes(pose):
            return None  # could alternately return empty set
        contacts = set()
        for res in pose:
            if res.has(self.nbr_atom):
                xyz = res.xyz(self.nbr_atom)
                nbrs = self.resi_map.neighbors(xyz)
                contacts.update(nbrs)
        return contacts


# below is only necessary to store stubs instead of resnums


def get_rosetta_stubs(pose, cen, aname1, aname2, aname3):
    centers, stubs = [], []
    for res in pose:
        if not (res.has(cen) and res.has(aname1) and res.has(aname2) and
                res.has(aname3)):
            continue
        centers.append(res.xyz(cen))
        stub = Stub(res.xyz(aname1), res.xyz(aname2), res.xyz(aname3))
        stubs.append((res.seqpos(), stub))
    return centers, stubs


class PoseStubIndex(object):
    def __init__(self, pose, clash_dis=3.2, nbr_dis=7.0, nbr_atom='CA'):
        bb_atoms = rcl.atoms(pose, 'bbheavy')
        cens, stubs = get_rosetta_stubs(pose, 'CA', 'N', 'CA', 'C')
        self.clash_check = stripe_index_3d(clash_dis, bb_atoms)
        self.stub_map = stripe_index_3d(nbr_dis, cens, stubs)
        self.nbr_atom = nbr_atom

    def clashes(self, pose):
        for res in pose:
            for ia in range(1, res.nheavyatoms() + 1):
                if self.clash_check.neighbor_exists(res.xyz(ia)):
                    return True
        return False

    def contacting_stubs(self, pose):
        if self.clashes(pose):
            return None  # could alternately return empty set
        contacts = set()
        for res in pose:
            if res.has(self.nbr_atom):
                nbrs = self.stub_map.neighbors(res.xyz(self.nbr_atom))
                contacts.update(nbrs)
        return contacts
