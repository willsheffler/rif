from __future__ import print_function
from rif.util import rcl
from rif.util.rcl import pyrosetta, rosetta
from rif.chem.biochem import aa_name3s
import os
import numpy as np
import pandas as pd


def print_full(x):
    pd.set_option('display.max_rows', len(x))
    print(x)
    pd.reset_option('display.max_rows')


def localpath(pth):
    assert not pth.startswith('/')
    return os.path.join(os.path.dirname(__file__), pth)


_cached_richardson_rotamer_space = None


def get_rotamer_space():
    """B is disulfide, x4 -> x2other"""
    global _cached_richardson_rotamer_space
    if _cached_richardson_rotamer_space is None:
        r = pd.read_csv(localpath('richardson.csv'), header=0, nrows=999, index_col=0)
        # print(r.columns)
        assert all(x > 0 for x in r.iloc[:, 0])
        for i in range(1, 5):
            r.iloc[:, i] = [float(x) / 100.0 for x in r.iloc[:, i]]
        assert all(0 <= r.f) and all(r.f <= 1.0)
        assert all(0 <= r.fa) and all(r.fa <= 1.0)
        assert all(0 <= r.fb) and all(r.fb <= 1.0)
        assert all(0 <= r.fo) and all(r.fo <= 1.0)
        # print_full(r.x1)
        assert sum(r.x1.isnull()) == 10
        assert all(-180.0 < r.x1.dropna()) and all(r.x1.dropna() <= 180.0)
        # print_full(r.x2)
        assert sum(r.x2.isnull()) == 15
        assert all(-180.0 < r.x2.dropna()) and all(r.x2.dropna() <= 180.0)
        # print_full(r.x3)
        assert sum(r.x3.isnull()) == 62
        assert all(-180.0 < r.x3.dropna()) and all(r.x3.dropna() <= 180.0)
        # print_full(r.x4)
        assert sum(r.x4.isnull()) == 92
        assert all(-180.0 < r.x4.dropna()) and all(r.x4.dropna() <= 180.0)
        assert sum(r.x1r.isnull()) == 10
        assert sum(r.x2r.isnull()) == 15
        # print_full(r.x3r)
        assert sum(r.x3r.isnull()) == 62
        # print_full(r.x4r)
        assert sum(r.x4r.isnull()) == 92
        for i in "1234":
            r['lb' + i] = np.NaN
            r['ub' + i] = np.NaN
            idx = r['x' + i + 'r'].notnull()
            r.loc[idx, 'lb' + i] = [float(x.split()[0]) for x in r['x' + i + 'r'][idx]]
            r.loc[idx, 'ub' + i] = [float(x.split()[2]) for x in r['x' + i + 'r'][idx]]
        # print_full(r.loc[:, 'lb1 ub1 ub2 ub2 lb3 ub3 lb4 ub4'.split()])
        # print_full(r.loc[:, 'x1w x2w x3w x4w'.split()])
        _cached_richardson_rotamer_space = r
    return _cached_richardson_rotamer_space


def sample_rotamer_space(rspace, resl=[10, 10, 10, 10]):
    """rotamer samples"""
    # for


def get_rotamer_index(rspace):
    """extract AA structural info via pyrosetta"""
    # print('get_rotamer_coords')
    rcl.init_check()
    chem_manager = rosetta.core.chemical.ChemicalManager
    rts = chem_manager.get_instance().residue_type_set("fa_standard")
    # print(rts)
    rfactory = rosetta.core.conformation.ResidueFactory
    for rname in aa_name3s:
        res = rfactory.create_residue(rts.name_map(rname))
        print(rname, res.nheavyatoms(), res.natoms())
    raise NotImplementedError
    return 1
