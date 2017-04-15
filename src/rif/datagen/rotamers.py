from __future__ import print_function
from rif.util import rcl
from rif.util.rcl import pyrosetta, rosetta
from rif.chem.biochem import aa_name3s

import pandas as pd


def get_rotamer_index():
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


def get_polar_heads():
    """return xarray.DataSets for atoms and raypairs"""
    atoms = list()
    rayxforms = list()
    raise NotImplementedError


richardson_rots = pd.DataFrame(
    [
        ('asn', 'p10', 60, 62, 0, 8, 20, -90, 0),
        ('asn', 'p30', 100, 62, 60, 6, 20, 0, 90),
        ('asn', 't20', 100, -174, -30, 5, 20, -120, 0),
        ('asn', 't30', 228, -177, 30, 14, 15, 0, 80),
        ('asn', 'm20', 580, -65, -20, 10, 20, -60, 10),
        ('asn', 'm80', 118, -65, -75, 10, 10, -100, -60),
        ('asn', 'm120', 58, -65, 120, 10, 18, 60, 100),
    ], columns='aa rot count chi1 chi2 sd1 sd2 lb2 ub2'.split()
)
