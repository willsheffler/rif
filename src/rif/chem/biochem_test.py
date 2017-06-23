import pytest

import rif
import rif.chem.biochem as bc
import rif.rcl as rcl


@pytest.mark.skipif('not rcl.HAVE_PYROSETTA')
def test_rif_atypes():
    rcl.init_check()
    rosetta_std_atypes = set()
    for resn in bc.aa_name3s:
        res = rcl.make_res(resn)
        for ia in range(res.natoms()):
            rosetta_std_atypes.add(str(res.atom_type(ia + 1).name()))
    for atype in rosetta_std_atypes:
        if atype == 'VIRT':
            continue
        assert atype in bc.rif_atype_names
    assert len(bc.aa_name1s) is len(bc.aa_name3s)


def test_chi_levers():
    chi_levers = bc.compute_chi_levers()
