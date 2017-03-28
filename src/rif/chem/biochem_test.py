import pytest

import rif
import rif.chem.biochem as biochem
import rif.util.rcl as rcl


@pytest.mark.skipif('not rcl.HAVE_PYROSETTA')
def test_rif_atypes():
    rcl.init_check()
    rosetta_std_atypes = set()
    for resn in biochem.aa_name3s:
        res = rcl.make_res(resn)
        for ia in range(res.natoms()):
            rosetta_std_atypes.add(str(res.atom_type(ia + 1).name()))
    for atype in rosetta_std_atypes:
        if atype == 'VIRT':
            continue
        assert atype in rif.chem.biochem.rif_atype_names
