from __future__ import print_function
from rif.util import rcl
from rif.util.rcl import pyrosetta, rosetta
from rif.chem.biochem import aa_name3s

import numpy as np
import pandas as pd
import xarray as xr


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
