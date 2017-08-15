import pytest
from rif import rcl
if rcl.HAVE_PYROSETTA:
    from rif.apps.examples.db_stub_index import PoseStubIndex, PoseResnumIndex


@pytest.mark.skipif('not rcl.HAVE_PYROSETTA')
def test_db_stub_index_1coi_clash(small_trimer_A, small_trimer_B):
    # small clash_dis 2.9 necessary for this test case
    cc = PoseStubIndex(small_trimer_A, clash_dis=2.9, nbr_atom='O')
    assert cc.contacting_stubs(small_trimer_A) is None  # clashes with self


@pytest.mark.skipif('not rcl.HAVE_PYROSETTA')
def test_db_stub_index_1coi_7(small_trimer_A, small_trimer_B):
    cc = PoseStubIndex(small_trimer_A, clash_dis=2.9,
                       nbr_dis=7.0, nbr_atom='O')
    assert len(cc.contacting_stubs(small_trimer_B)) == 0


@pytest.mark.skipif('not rcl.HAVE_PYROSETTA')
def test_db_stub_index_1coi_765(small_trimer_A, small_trimer_B):
    cc = PoseStubIndex(small_trimer_A, clash_dis=2.9,
                       nbr_dis=7.65, nbr_atom='O')
    assert len(cc.contacting_stubs(small_trimer_B)) == 5


@pytest.mark.skipif('not rcl.HAVE_PYROSETTA')
def test_db_stub_index_1coi_10(small_trimer_A, small_trimer_B):
    cc = PoseStubIndex(small_trimer_A, clash_dis=2.9,
                       nbr_dis=10.0, nbr_atom='O')
    assert len(cc.contacting_stubs(small_trimer_B)) == 13


@pytest.mark.skipif('not rcl.HAVE_PYROSETTA')
def test_db_stub_index_1coi_20(small_trimer_A, small_trimer_B):
    cc = PoseStubIndex(small_trimer_A, clash_dis=2.9,
                                         nbr_dis=20.0, nbr_atom='O')
    assert len(cc.contacting_stubs(small_trimer_B)) == 29


@pytest.mark.skipif('not rcl.HAVE_PYROSETTA')
def test_db_resi_index_1coi_clash(small_trimer_A, small_trimer_B):
    # small clash_dis 2.9 necessary for this test case
    cc = PoseResnumIndex(small_trimer_A, clash_dis=2.9, nbr_atom='O')
    assert cc.contacting_resnums(small_trimer_A) is None  # clashes with self


@pytest.mark.skipif('not rcl.HAVE_PYROSETTA')
def test_db_resi_index_1coi_7(small_trimer_A, small_trimer_B):
    cc = PoseResnumIndex(small_trimer_A, clash_dis=2.9,
                         nbr_dis=7.0, nbr_atom='O')
    assert len(cc.contacting_resnums(small_trimer_B)) == 0


@pytest.mark.skipif('not rcl.HAVE_PYROSETTA')
def test_db_resi_index_1coi_765(small_trimer_A, small_trimer_B):
    cc = PoseResnumIndex(small_trimer_A, clash_dis=2.9,
                         nbr_dis=7.65, nbr_atom='O')
    assert cc.contacting_resnums(small_trimer_B) == {5, 8, 12, 15, 19}


@pytest.mark.skipif('not rcl.HAVE_PYROSETTA')
def test_db_resi_index_1coi_10(small_trimer_A, small_trimer_B):
    cc = PoseResnumIndex(small_trimer_A, clash_dis=2.9,
                         nbr_dis=10.0, nbr_atom='O')
    contacts = cc.contacting_resnums(small_trimer_B)
    assert contacts == {1, 2, 5, 8, 9, 12, 15, 16, 19, 20, 22, 23, 26}


@pytest.mark.skipif('not rcl.HAVE_PYROSETTA')
def test_db_resi_index_1coi_20(small_trimer_A, small_trimer_B):
    cc = PoseResnumIndex(small_trimer_A, clash_dis=2.9,
                         nbr_dis=20.0, nbr_atom='O')
    assert len(cc.contacting_resnums(small_trimer_B)) == 29
