
from builtins import *

import random
from math import sqrt

from collections import defaultdict
import sys
import multiprocessing
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import pandas as pd

import rif
from rif import rcl
from rif.rcl import AtomID


def make_atype_charges():
    rcl.init_check()
    avg_iatype_charge = dict()
    atomic_charges = defaultdict(list)
    for resn in rif.biochem.aa_name3s:
        res = rcl.make_res(resn)
        for ia in range(res.natoms()):
            iatype = res.atom_type_index(ia + 1)
            # if iatype > 100:
            #     print iatype, resn, res.atom_name(ia+1)
            charge = res.atomic_charge(ia + 1)
            atomic_charges[iatype].append(charge)
    for i, charges in list(atomic_charges.items()):
        # print i, rcl.ats()[i].name(), np.mean(charges), np.std(charges)
        avg_iatype_charge[rcl.ats()[i].name()] = np.mean(charges)
    return avg_iatype_charge


def patch_database_with_single_atom_params_files(mock_atom='MOCK'):
    mock_atom = mock_atom.upper()
    new_atom_types_txt = """
## mock atom which *should* never have any score, but otherwise
## act as a normal heavy atom (will allow attached Hs to be scored)
MOCK    C    0.0000    0.0000    0.0000    3.5000    0.0000
"""
    new_extras = {
        'soft_rep_params': 'MOCK ' + '0 ' * 1 + ' ## RIF\n',
        'gen_born_params': 'MOCK ' + '0 ' * 2 + ' ## RIF\n',
        'gen_kirkwood_params': 'MOCK ' + '0 ' * 1 + ' ## RIF\n',
        'sasa_radii_legacy': 'MOCK ' + '0 ' * 1 + ' ## RIF\n',
        'NACCESS_sasa_radii': 'MOCK ' + '0 ' * 1 + ' ## RIF\n',
        'reduce_sasa_radii': 'MOCK ' + '0 ' * 1 + ' ## RIF\n',
        'memb_fa_params': 'MOCK ' + '0 ' * 4 + ' ## RIF\n',
        'std_charges': 'MOCK ' + '0 ' * 1 + ' ## RIF\n',
        'atom_orbital_hybridization': 'MOCK ' + '0 ' * 3 + ' ## RIF\n',
        'facts_born_params.epm1': 'MOCK ' + '0 ' * 17 + ' ## RIF\n',
    }

    singe_atom_res_template = """NAME {name}
IO_STRING  {name} Z
TYPE LIGAND
AA UNK
ATOM  A1   {mock} X    0.0
ATOM  A2   VIRT  X    0.0
ATOM  A3   {name} {mmname} {charge}
BOND A1 A2
BOND A1 A3

ICOOR_INTERNAL A1 0.0  0.0 0.0 A1 A2 A3
ICOOR_INTERNAL A2 0.0  0.0 1.0 A1 A2 A3
ICOOR_INTERNAL A3 0.0 90.0 1.0 A1 A2 A3

PROPERTIES
NBR_ATOM A3
NBR_RADIUS 8.0
FIRST_SIDECHAIN_ATOM ALL
"""
    up_to_date = True
    atype_charges = make_atype_charges()
    dbpath = 'chemical/residue_type_sets/fa_standard/residue_types/single_atoms/'

    for atype in rif.biochem.rif_atype_names:
        content = singe_atom_res_template.format(
            name=atype, charge=atype_charges[atype], mmname='X', mock=mock_atom)
        db_changed = rif.rcl.add_files_to_rosettaDB(dbpath + atype + '.params', content)
        up_to_date = up_to_date and not db_changed

    params_file = 'residue_types/single_atoms/%s.params'
    params_files = [params_file % x for x in rif.biochem.rif_atype_names]
    db_changed = rif.rcl.add_residue_types_to_rosettaDB(params_files, '## single atoms')
    up_to_date = up_to_date and not db_changed

    db_changed = rcl.add_atom_types_to_rosettaDB(new_atom_types_txt, new_extras)
    up_to_date = up_to_date and not db_changed

    if not up_to_date:
        print('!' * 80)
        print('!' * 23 + ' RosettaDB UPDATED, PLEASE RE-RUN ' + '!' * 23)
        print('!' * 80)
        sys.exit()


def identity(x): return x
def square(x): return x * x


class EtableSpec(object):
    """metadata for an Etable"""
    def __init__(self, mn=0.0, mx=6.0, n=32, step=None,
                 map=(identity, identity), label='beta_nov15'):
        if not n:
            n = int((self.mx - self.mn) / self.step)
        self.map = map
        self.mn = map[0](mn)
        self.mx = map[0](mx)
        self.n = n
        self.step = (self.mx - self.mn) / n
        self.label = label
        dists = []
        for i in range(self.n):
            dists.append(map[1](float(i) * self.step + self.mn))
        self.dists = np.array(dists)


def zero_residues_on_atom(pose, atomno):
    for ir in range(1, pose.n_residue() + 1):
        pos0 = pose.residue(ir).xyz(atomno)
        for ia in range(1, pose.residue_type(ir).natoms() + 1):
            pose.set_xyz(AtomID(ia, ir), pose.xyz(AtomID(ia, ir)) - pos0)


def run_make_eframe_for_atoms(args):
    return make_eframe_for_atoms(*args)


def make_eframe_for_atoms(iat1, iat2, atype1, atype2, spec):
    sfunc = rcl.create_score_function(spec.label)
    stypes = sfunc.get_nonzero_weighted_scoretypes()
    stypes = [stypes[i] for i in range(1, len(stypes) + 1)]
    iref = 3
    pose = rcl.Pose()
    pose.append_residue_by_jump(rcl.make_res(atype1), 0)
    pose.append_residue_by_jump(rcl.make_res(atype2), 1)
    zero_residues_on_atom(pose, iref)
    assert pose.residue(1).xyz(iref).length() < 0.001
    assert pose.residue(2).xyz(iref).length() < 0.001
    assert pose.num_jump() == 1
    orig_jump_trans = rcl.xyzVector(pose.jump(1).get_translation())  # copy
    trans = pose.residue(2).xyz(iref) - pose.residue(1).xyz(iref)
    if trans.length() < 0.001:
        trans = rcl.xyzVector(1, 0, 0)
    trans.normalize()
    eframe = list()
    w = sfunc.weights()
    for i, d in enumerate(spec.dists):
        j = pose.jump(1)
        d_trans = rcl.xyzVector(d * trans.x, d * trans.y, d * trans.z)
        j.set_translation(orig_jump_trans + d_trans)
        pose.set_jump(1, j)
        dpose = pose.residue(1).xyz(iref).distance(pose.residue(2).xyz(iref))
        assert abs(dpose - d) < 0.001
        score = sfunc(pose)
        emap = pose.energies().total_energies()
        assert abs(emap.dot(w) - score) < 0.001
        getname = rcl.core.scoring.name_from_score_type
        scores = {str(getname(st)): emap.get(st) * w.get(st) for st in stypes}
        scores['score'] = score
        scores['atype1'] = atype1
        scores['atype2'] = atype2
        scores['dist'] = d
        thisframe = pd.DataFrame(scores, index=[i])
        eframe.append(thisframe)
    eframe = pd.concat(eframe)
    # eframe = eframe.set_index(['atype1', 'atype2', 'dist'], drop=True)
    # eframe = eframe.astype(np.float32)
    sys.stdout.write('.')
    if 0 == random.randint(0, 50):
        sys.stdout.write('\n')
    sys.stdout.flush()
    return eframe


def make_eframes(spec, multiprocess=False,
                 atypes1=rif.biochem.rif_atype_names,
                 atypes2=rif.biochem.rif_atype_names):
    # type: (EtableSpec, bool, list[str], list[str]) -> pd.DataFrame
    jobs = list()
    for iat1, atype1 in enumerate(atypes1):
        for iat2, atype2 in enumerate(atypes2):
            if atypes1 == atypes2 and iat1 > iat2:
                continue
            jobs.append((iat1, iat2, atype1, atype2, spec))
    if multiprocess:
        ncpu = multiprocessing.cpu_count() // 2
        with ProcessPoolExecutor(ncpu) as executor:
            eframes = executor.map(run_make_eframe_for_atoms, jobs)
    else:
        eframes = list(map(run_make_eframe_for_atoms, jobs))
    eframe = pd.concat(eframes)
    return eframe


if __name__ == '__main__':
    rcl.init_check()
    mock_atom = 'mock'
    patch_database_with_single_atom_params_files(mock_atom=mock_atom)
    spec = EtableSpec(mn=0.0, mx=6, n=64, label='beta_nov15', map=(square,sqrt))
    # print(spec.dists)
    # sys.exit(0)
    sfunc = rcl.create_score_function(spec.label)
    eframe = make_eframes(spec, multiprocess=True)
    # atypes1='CH3 Hapo'.split(),
    # atypes2='CH3 OOC Nbb Hpol Hapo Haro HNbb'.split())
    print('made eframe', eframe.shape)
    path = rif.resources.locate_resource_file('rosetta_score_data.h5')
    h5path = '/eframes/' + mock_atom + '_dsq' + '/beta_nov15/eframe'
    print('save to', path, ':', h5path)
    store = pd.HDFStore(path)
    store[h5path] = eframe
    store.close()
