import pytest

from rif.util import rcl


@pytest.mark.skipif('not rcl.HAVE_PYROSETTA')
def test_import_pyrosetta():
    import pyrosetta
    import rosetta


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
