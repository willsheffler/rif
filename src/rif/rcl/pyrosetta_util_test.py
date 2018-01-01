import pytest
from types import ModuleType
from mock import MagicMock
from rif import rcl


def test_import_pyrosetta():
    if rcl.HAVE_PYROSETTA:
        assert isinstance(rcl.pyrosetta, ModuleType)
        assert isinstance(rcl.rosetta, ModuleType)


@pytest.mark.skipif('not rcl.HAVE_PYROSETTA')
def test_pyrosetta_init():
    rcl.init_check()
    rcl.init_check()  # second call ignored
    with pytest.raises(rcl.ReInitError) as e:
        rcl.init_check('-score:weighs foo')  # change options fails
    assert 'previous' in str(e)
    assert 'thiscall' in str(e)


@pytest.mark.skipif('not rcl.HAVE_PYROSETTA')
def test_make_ideal_res():
    rcl.init_check()
    res = rcl.make_res('ALA')
    assert isinstance(res, rcl.core.conformation.Residue)
    assert res.name3() == 'ALA'


@pytest.mark.skipif('not rcl.HAVE_PYROSETTA')
def test_atomic_charge_deviation():
    pass


# @pytest.mark.skipif('not rcl.HAVE_PYROSETTA')
# def test_generate_canonical_residue():
#     rcl.init_check()
#     k = rcl.generate_canonical_residue('LYS')
#     assert k.name1() == 'K'
#     assert k.nheavyatoms() == 9


# @pytest.mark.skipif('not rcl.HAVE_PYROSETTA')
# def test_generate_canonical_rotamer_residues():
#     rcl.init_check()
#     rots = rcl.generate_canonical_rotamer_residues('ASN')
#     assert len(rots) in (267, 399)  # WTF Rosetta!
#     # assert len(rots) == 399
#     # for i, rot in enumerate(rots):
#     # p = rcl.Pose()
#     # p.append_residue_by_jump(rot, 0)
#     # print(i)
#     # p.dump_pdb('asn%i.pdb' % i)
#     # print(len(rots))
#     # print(rcl.rosetta.options)
