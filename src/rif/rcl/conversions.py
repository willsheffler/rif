from . import pyrosetta_util
import rif
from rif.sele import parse_atom_names, parse_ray_atom_names
from rif.chem.biochem import rif_atype_names
from rif import Atom
from rif.geom import Ray
import numpy as np
from rif import rcl


_RIF_ATYPE_MAP = None

_RIF_ATOM_NAME_MAP = dict()
for aname in rif_atype_names:
    _RIF_ATOM_NAME_MAP[aname] = aname
del aname
_RIF_ATOM_NAME_MAP['CH0'] = 'aroC'
_RIF_ATOM_NAME_MAP['HS'] = 'Hpol'
_RIF_ATOM_NAME_MAP['NtrR'] = 'Narg'
_RIF_ATOM_NAME_MAP['SH1'] = 'S'


def rosetta2rif_atom_name(rosetta_atom_name):
    return _RIF_ATOM_NAME_MAP[rosetta_atom_name]


def get_rif_atype_map():
    global _RIF_ATYPE_MAP
    if None is _RIF_ATYPE_MAP:
        pyrosetta_util.init_check(strict=False)
        for aname in _RIF_ATOM_NAME_MAP.values():
            assert aname in rif_atype_names
        ats = pyrosetta_util.ats()
        atypemap = np.zeros(ats.n_atomtypes() + 1, dtype='i4')
        atypemap[:] = 255
        for aname in _RIF_ATOM_NAME_MAP.keys():
            rif_atype_index = rif_atype_names.index(
                rosetta2rif_atom_name(aname))
            atypemap[ats.atom_type_index(aname)] = rif_atype_index
        _RIF_ATYPE_MAP = atypemap
    return _RIF_ATYPE_MAP


def rif_atype(rosetta_atype):
    return get_rif_atype_map()[rosetta_atype]


def rif_rtype(rosetta_aa):
    return int(rosetta_aa) - 1


def _get_atom_from_res(res, aname, verbose=False):
    if not res.has(aname):
        raise AttributeError(
            'pose res %i %s has no atom %s' % (res.seqpos(), res.name3(), aname))
    xyz = res.xyz(aname)
    atomno = res.atom_index(aname)
    if verbose:
        print('_get_atom_from_res:', res.name3(), atomno,
              res.atom_name(atomno),
              res.atom_type_index(atomno),
              rif_atype(res.atom_type_index(atomno)))
    atype = rif_atype(res.atom_type_index(atomno))
    restype = rif_rtype(res.aa())
    return (
        ([xyz.x, xyz.y, xyz.z],),
        atype,
        atomno - 1,  # 0-indexing!
        restype,
    )


def first_available_atom(res, anames):
    for aname in anames:
        if res.has(aname):
            return aname
    return None


def atoms_fixed_width(pose, sele, protein_only=True, **kwargs):
    anames_list = parse_atom_names(sele)
    atoms = list()
    for ir in range(1, pose.size() + 1):
        res = pose.residue(ir)
        if protein_only and not res.is_protein():
            continue
        for anames in anames_list:
            aname = first_available_atom(res, anames)
            atoms.append(_get_atom_from_res(res, aname, **kwargs))
    atoms = np.array(atoms, dtype=Atom)
    n = len(anames_list)
    if n > 1:
        assert len(atoms) % n == 0
        atoms = atoms.reshape(int(len(atoms) / n), n)
    return atoms


def atoms_matching(pose, predicate, **kwargs):
    atoms = list()
    for ir in range(1, pose.size() + 1):
        res = pose.residue(ir)
        for ia in range(1, res.natoms() + 1):
            aname = res.atom_name(ia)
            if predicate(res, aname):
                atoms.append(_get_atom_from_res(res, aname, **kwargs))
    atoms = np.array(atoms, dtype=Atom)
    return atoms


def pred_all(*args, **kwargs):
    return True


def pred_heavy(res, aname, **kwargs):
    return not res.atom_is_hydrogen(res.atom_index(aname))


def pred_bb(res, aname, **kwargs):
    return res.atom_is_backbone(res.atom_index(aname))


def pred_bbheavy(res, aname, **kwargs):
    return (res.atom_is_backbone(res.atom_index(aname)) and
            not res.atom_is_hydrogen(res.atom_index(aname)))


def atoms(pose, sele='HEAVY', **kwargs):
    'extract rif style atoms from rosetta pose'
    if sele.lower() == 'all':
        return atoms_matching(pose, pred_all, **kwargs)
    elif sele.lower() == 'heavy':
        return atoms_matching(pose, pred_heavy, **kwargs)
    elif sele.lower() == 'bb':
        return atoms_matching(pose, pred_bb, **kwargs)
    elif sele.lower() == 'bbheavy':
        return atoms_matching(pose, pred_bbheavy, **kwargs)
    else:
        return atoms_fixed_width(pose, sele=sele, **kwargs)


# why this instead of rays(atoms(from), atoms(to))?
def rays(pose, sele, shifts=None, protein_only=True, **kwargs):
    'extract rif style rays from rosetta pose'
    ray_list = parse_ray_atom_names(sele)
    if shifts is None:
        shifts = [0] * len(ray_list)
    rays = list()
    for ir in range(1, pose.size() + 1):
        rays_tmp = []
        # tmp = 0
        for anames, shift in zip(ray_list, shifts):
            if ir + shift > pose.size():
                break
            res = pose.residue(ir + shift)
            if protein_only and not res.is_protein():
                break
            aname0 = first_available_atom(res, anames[0])
            aname1 = first_available_atom(res, anames[1])
            if not aname0 or not aname1:
                break
            orig = res.xyz(aname0)
            dest = res.xyz(aname1)
            dirn = (dest - orig).normalized()
            ray = ((((orig.x, dirn.x),
                     (orig.y, dirn.y),
                     (orig.z, dirn.z),
                     (1, 0),),),)
            # ray = ((((0, 0),
            # (0, 0),
            # (0, 0),
            # (ir + shift, tmp),),),)
            # tmp += 1
            rays_tmp.append(ray)
        else:
            rays.extend(rays_tmp)
    rays = np.array(rays, dtype=Ray)
    n = len(ray_list)
    if n > 1:
        assert len(rays) % n == 0
        rays = rays.reshape(int(len(rays) / n), n)
    return rays


def bbstubs(pose, which_resi=None):
    'extract rif style stubs from rosetta pose'
    if which_resi is None:
        which_resi = list(range(1, pose.size() + 1))
    n_prot_res = 0
    for ir in which_resi:
        n_prot_res += pose.residue(ir).is_protein()
    rif_stubs = np.zeros(n_prot_res, dtype='(4,4)f')
    n_prot_res = 0
    for r in (pose.residue(i) for i in which_resi):
        if not r.is_protein():
            continue
        ros_stub = rcl.Stub(r.xyz('CA'), r.xyz('N'), r.xyz('CA'), r.xyz('C'))
        rif_stub = to_rif_stub(ros_stub)
        rif_stubs[n_prot_res, :, :] = rif_stub['raw']
        n_prot_res += 1
    return rif_stubs.reshape(n_prot_res * 16).view(rif.X3)


def to_rif_stub(rosstub):
    rifstub = np.zeros(1, dtype='4,4f')
    for i in range(3):
        rifstub[..., i, 3] = rosstub.v[i]
        for j in range(3):
            rifstub[..., i, j] = rosstub.M(i + 1, j + 1)
    rifstub[..., 3, 3] = 1.0
    return rifstub.reshape(16,).view(rif.X3)


def to_rosetta_stub(rifstub, i=0):
    if isinstance(rifstub, rif.X3):
        rifstub = np.array([(rifstub,)], dtype=rif.X3)
    if isinstance(i, int):
        i = (i,)
    if rifstub.shape[-2:] != (4, 4):
        rifstub = rifstub.view('4,4f')
    if rifstub.shape == (4, 4):
        rifstub = rifstub[np.newaxis, ...]
    rifstub = rifstub.astype('f4')
    rosstub = rcl.Stub()
    rosstub.M.xx = rifstub[i + (0, 0)]
    rosstub.M.xy = rifstub[i + (0, 1)]
    rosstub.M.xz = rifstub[i + (0, 2)]
    rosstub.M.yx = rifstub[i + (1, 0)]
    rosstub.M.yy = rifstub[i + (1, 1)]
    rosstub.M.yz = rifstub[i + (1, 2)]
    rosstub.M.zx = rifstub[i + (2, 0)]
    rosstub.M.zy = rifstub[i + (2, 1)]
    rosstub.M.zz = rifstub[i + (2, 2)]
    rosstub.v.x = rifstub[i + (0, 3)]
    rosstub.v.y = rifstub[i + (1, 3)]
    rosstub.v.z = rifstub[i + (2, 3)]
    return rosstub
