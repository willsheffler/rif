

"""
Rosetta Compatibility Layer

if pyrosetta is available, provides pyrosetta. check value of HAVE_PYROSETTA.
"""

import os

try:
    import pyrosetta
    import rosetta
    from rosetta import utility, numeric, basic, core
    from rosetta.core.pose import make_pose_from_sequence, Pose
    from rosetta.core.kinematics import FoldTree, MoveMap
    from rosetta.core.import_pose import pose_from_file
    from rosetta.core.io.pdb import dump_pdb
    from rosetta.core.id import AtomID
    from rosetta.core.scoring import ScoreFunction, get_score_function
    from rosetta.numeric import xyzVector_double_t as xyzVector

    create_score_function = rosetta.core.scoring.ScoreFunctionFactory.create_score_function
    ROSETTA_DB = pyrosetta._rosetta_database_from_env()
    HAVE_PYROSETTA = True
except ImportError as e:
    import mock
    HAVE_PYROSETTA = False
    pyrosetta = mock.MagicMock()
    rosetta = mock.MagicMock()




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


def init_check(options='-beta_nov15', fail_if_options_dont_match=True):
    global pyrosetta_init_options
    if pyrosetta_init_options is None:
        pyrosetta.init(options=options)
        pyrosetta_init_options = options
    elif options != pyrosetta_init_options:
        if fail_if_options_dont_match:
            raise ReInitError(
                'attempt to init rosetta with different options: previous: {}, thiscall: {}'
                .format(pyrosetta_init_options, options))


def ats():
    return core.chemical.ChemicalManager.get_instance().atom_type_set('fa_standard')


def rts():
    return core.chemical.ChemicalManager.get_instance().residue_type_set('fa_standard')


def get_rtype(res_type_name):
    return rts().name_map(res_type_name)


def make_res(res_type_name):
    return core.conformation.ResidueFactory.create_residue(get_rtype(res_type_name))


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
            f.write("\n\n" + comment + '\n')
            for newline in update_needed:
                f.write(newline + '\n')
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
