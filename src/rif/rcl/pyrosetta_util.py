

"""
nRosetta Compatibility Layer

if pyrosetta is available, provides pyrosetta. check value of HAVE_PYROSETTA.
"""

import os

try:
    import pyrosetta
    from pyrosetta import _rosetta_database_from_env
    try:
        from pyrosetta import rosetta
        from pyrosetta.rosetta import utility, numeric, basic, core
        from pyrosetta.rosetta.core.conformation import Residue
        from pyrosetta.rosetta.core.pose import make_pose_from_sequence, Pose
        from pyrosetta.rosetta.core.kinematics import FoldTree, MoveMap, Stub
        from pyrosetta.rosetta.core.import_pose import pose_from_file
        from pyrosetta.rosetta.core.io.pdb import dump_pdb
        from pyrosetta.rosetta.core.id import AtomID
        from pyrosetta.rosetta.core.scoring import ScoreFunction, get_score_function
        from pyrosetta.rosetta.numeric import xyzVector_double_t as xyzVector
        from pyrosetta.rosetta.numeric import xyzVector_double_t as xyzMatrix
        import pyrosetta.rosetta.basic
        import pyrosetta.rosetta.core
        import pyrosetta.rosetta.core.pack.rotamer_set
        from pyrosetta.rosetta.protocols.protein_interface_design.movers import TryRotamers
    except ImportError:
        import rosetta
        from rosetta import utility, numeric, basic, core
        from rosetta.core.conformation import Residue
        from rosetta.core.pose import make_pose_from_sequence, Pose
        from rosetta.core.kinematics import FoldTree, MoveMap, Stub
        from rosetta.core.import_pose import pose_from_file
        from rosetta.core.io.pdb import dump_pdb
        from rosetta.core.id import AtomID
        from rosetta.core.scoring import ScoreFunction, get_score_function
        from rosetta.numeric import xyzVector_double_t as xyzVector
        from rosetta.numeric import xyzVector_double_t as xyzMatrix
        import rosetta.basic
        import rosetta.core
        import rosetta.core.pack.rotamer_set
        from rosetta.protocols.protein_interface_design.movers import TryRotamers

    create_score_function = rosetta.core.scoring.ScoreFunctionFactory.create_score_function
    # for x in dir(pyrosetta):
    # print('pyrosetta:', x)/
    ROSETTA_DB = _rosetta_database_from_env()
    HAVE_PYROSETTA = True
except ImportError as e:
    # import mock
    HAVE_PYROSETTA = False
    # pyrosetta = mock.MagicMock()
    # rosetta = mock.MagicMock()

from .conversions import to_rosetta_stub, bbstubs
import numpy as np


class Error(Exception):
    """Rosetta Compatibility Exception"""


class ReInitError(Error):
    """Tried to reinitialize rosetta with different options"""


class AlreadyInitError(Error):
    """Rosetta already initialized but should not be"""


class RosettaDBError(Error):
    """Rosetta Database Error"""


pyrosetta_init_options = None


def is_initialized():
    return pyrosetta_init_options is not None


def init_check(options=None, strict=True):
    if options is None:
        strict = False
        options = '-corrections:beta_nov16 -mute all'
    global pyrosetta_init_options
    if pyrosetta_init_options is None:
        pyrosetta.init(options=options)
        pyrosetta_init_options = options
    elif options != pyrosetta_init_options:
        if strict:
            raise ReInitError(
                'attempt to init rosetta with different options: previous: {}, thiscall: {}'
                .format(pyrosetta_init_options, options))


def ats(kind='fa_standard'):
    init_check()
    return core.chemical.ChemicalManager.get_instance().atom_type_set(kind)


def rts(kind='fa_standard'):
    init_check()
    return core.chemical.ChemicalManager.get_instance().residue_type_set(kind)


def get_rtype(res_type_name):
    return rts().name_map(res_type_name)


def make_res(res_type_name):
    return core.conformation.ResidueFactory.create_residue(
        get_rtype(res_type_name))


def add_files_to_rosettaDB(fname, contents):
    """returns true if update was done"""
    mypath = ROSETTA_DB + '/' + fname
    if not os.path.exists(os.path.dirname(mypath)):
        os.mkdir(os.path.dirname(mypath))
    update_needed = False
    if os.path.exists(mypath):
        with open(mypath, 'r') as f:
            if f.read() != contents:
                update_needed = True
                # raise RosettaDBError("patched rosetta database " \
                #                      "file exists and dosen't match: " + mypath)
                print('updating rosettaDB:', fname)
    if update_needed:
        with open(mypath, 'w') as f:
            f.write(contents)
        print('added to rosettaDB:', mypath)
    return update_needed


def add_residue_types_to_rosettaDB(newfiles, comment="## RIF added restypes",
                                   typeset='fa_standard'):
    """return true if update was done"""
    assert newfiles
    mypath = (ROSETTA_DB + '/chem/residue_type_sets/' +
              typeset + '/residue_types.txt')
    update_needed = list()
    with open(mypath, 'r') as f:
        existing_files = set(x.strip() for x in f)
        for newfile in newfiles:
            if newfile not in existing_files:
                update_needed.append(newfile)
    if update_needed:
        with open(mypath, 'a') as f:
            f.write(os.linesep * 2 + comment + os.linesep)
            for newline in update_needed:
                f.write(newline + os.linesep)
        print('modified rosettaDB:', mypath)
    else:
        print('nochange rosettaDB:', mypath)
    return bool(update_needed)


def update_roesettaDB_file(fname, newtext):
    if not os.path.exists(fname):
        raise RosettaDBError('rosettaDB file doesnt exist: ' + fname)
    with open(fname, 'r') as f:
        update_needed = newtext not in f.read()
    if update_needed:
        with open(fname, 'a') as f:
            f.write(newtext)
        print('modified rosettaDB:', fname)
    else:
        print('nochange rosettaDB:', fname)
    return update_needed


def add_atom_types_to_rosettaDB(newprops, extradict, typeset='fa_standard'):
    any_updates = False
    atomprops = (ROSETTA_DB + '/chem/atom_type_sets/' +
                 typeset + '/atom_properties.txt')
    updated = update_roesettaDB_file(atomprops, newprops)
    any_updates = any_updates or updated
    for extraname, extraval in list(extradict.items()):
        mypath = (ROSETTA_DB + '/chem/atom_type_sets/' +
                  typeset + '/extras/' + extraname + '.txt')
        updated = update_roesettaDB_file(mypath, extraval)
        any_updates = any_updates or updated
    return any_updates


def generate_canonical_residue(residue_name3):
    work_pose = Pose()
    make_pose_from_sequence(
        work_pose, "X[%s]" % residue_name3, "fa_standard", auto_termini=False)
    work_residue = rosetta.core.conformation.Residue(work_pose.residue(1))

    ca_loc = work_residue.xyz("CA")

    for a in range(work_residue.natoms()):
        work_residue.set_xyz(a + 1, work_residue.xyz(a + 1) - ca_loc)

    return work_residue


def generate_canonical_rotamer_residues_phipsi(residue_name3, target_phi_psi):
    raise NotImlementedError
    canonical_residue = generate_canonical_residue(residue_name3)
    test_sequence = "AAX[%s]AA" % residue_name3
    target_phi, target_psi = target_phi_psi
    sf = get_score_function()
    tryrot = TryRotamers(
        resnum=3,
        scorefxn=sf,
        explosion=0,
        jump_num=0,
        clash_check=True,
        solo_res=False,
        include_current=False)
    test_pose = Pose()
    make_pose_from_sequence(test_pose, test_sequence, "fa_standard")
    for i in range(1, test_pose.size() + 1):
        test_pose.set_psi(i, target_psi)
        test_pose.set_phi(i, target_phi)
        test_pose.set_omega(i, 180)
    tryrot.setup_rotamer_set(test_pose)
    rotamer_set = tryrot.rotamer_set()
    # print('rotamer_set.num_rotamers()',
    # residue_name3, target_phi_psi, rotamer_set.num_rotamers())
    rotamers = [rotamer_set.rotamer(i).clone()
                for i in range(1, rotamer_set.num_rotamers() + 1)]
    for r in rotamers:
        r.orient_onto_residue(canonical_residue)
        r.seqpos(1)
    return rotamers


def generate_canonical_rotamer_residues(residue_name3):
    canonical_phi_psi = {
        "helical": (-66.0, -42.0),
        "sheet": (-126.0, 124.0),
    }
    result = list()
    for phipsi in canonical_phi_psi.values():
        result.extend(
            generate_canonical_rotamer_residues_phipsi(residue_name3, phipsi))
        # print(phipsi, len(result))
    return result


def worst_CN_connect(p):
    for ir in range(1, len(p)):
        worst = 0
        if (p.residue(ir).is_protein() and
                p.residue(ir + 1).is_protein() and not (
                rosetta.core.pose.is_upper_terminus(p, ir) or
                rosetta.core.pose.is_lower_terminus(p, ir + 1))):
            dist = p.residue(ir).xyz('C').distance(p.residue(ir + 1).xyz('N'))
            worst = max(abs(dist - 1.32), worst)
    return worst


def pose_bounds(pose, lb, ub):
    if ub < 0: ub = len(pose) + 1 + ub
    if lb < 1 or ub > len(pose):
        raise ValueError('lb/ub ' + str(lb) + '/' + str(ub) +
                         ' out of bounds for pose with len '
                         + str(len(pose)))
    return lb, ub


def xform_pose(xform, pose, lb=1, ub=-1):
    lb, ub = pose_bounds(pose, lb, ub)
    if xform.shape != (4, 4):
        raise ValueError(
            'invalid xform, must be 4x4 homogeneous matrix, shape is: '
            + str(xform.shape))
    xform = to_rosetta_stub(xform)
    rosetta.protocols.sic_dock.xform_pose(pose, xform, lb, ub)


def concatenate_pose(p, to_append, lb=1, ub=-1, out=None):
    "put poses together aligning on end of p"
    lb, ub = pose_bounds(to_append, lb, ub)
    if not p.residue(len(p)).is_protein():
        raise ValueError('last res of p is not protein')
    if not to_append.residue(lb).is_protein():
        raise ValueError('lb res of to_append is not protein')
    if out is not p:
        out = p.clone()
    stub1 = bbstubs(p, [len(p)])['raw'][0]
    stub2 = bbstubs(to_append, [lb])['raw'][0]
    xalign = stub1 @ np.linalg.inv(stub2)
    rosetta.core.pose.remove_upper_terminus_type_from_pose_residue(
        out, len(out))
    outlen = len(out)
    rosetta.core.pose.append_subpose_to_pose(
        out, to_append, lb + 1, ub, new_chain=0)
    xform_pose(xalign, out, outlen + 1)
    # todo: fix NH/CO
    return out
