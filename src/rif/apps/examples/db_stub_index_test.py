import pytest
from rif import rcl
if rcl.HAVE_PYROSETTA:
    from rif.apps.examples.db_stub_index import PoseStubIndex, PoseResnumIndex


@pytest.mark.skipif('not rcl.HAVE_PYROSETTA')
def test_db_stub_index_1coi_clash(trimerA_pose, trimerB_pose):
    # small clash_dis 2.9 necessary for this test case
    cc = PoseStubIndex(trimerA_pose, clash_dis=2.9, nbr_atom='O')
    assert cc.contacting_stubs(trimerA_pose) is None  # clashes with self


@pytest.mark.skipif('not rcl.HAVE_PYROSETTA')
def test_db_stub_index_1coi_7(trimerA_pose, trimerB_pose):
    cc = PoseStubIndex(trimerA_pose, clash_dis=2.9,
                       nbr_dis=7.0, nbr_atom='O')
    assert len(cc.contacting_stubs(trimerB_pose)) == 0


@pytest.mark.skipif('not rcl.HAVE_PYROSETTA')
def test_db_stub_index_1coi_765(trimerA_pose, trimerB_pose):
    cc = PoseStubIndex(trimerA_pose, clash_dis=2.9,
                       nbr_dis=7.65, nbr_atom='O')
    assert len(cc.contacting_stubs(trimerB_pose)) == 5


@pytest.mark.skipif('not rcl.HAVE_PYROSETTA')
def test_db_stub_index_1coi_10(trimerA_pose, trimerB_pose):
    cc = PoseStubIndex(trimerA_pose, clash_dis=2.9,
                       nbr_dis=10.0, nbr_atom='O')
    assert len(cc.contacting_stubs(trimerB_pose)) == 13


@pytest.mark.skipif('not rcl.HAVE_PYROSETTA')
def test_db_stub_index_1coi_20(trimerA_pose, trimerB_pose):
    cc = PoseStubIndex(trimerA_pose, clash_dis=2.9,
                       nbr_dis=20.0, nbr_atom='O')
    assert len(cc.contacting_stubs(trimerB_pose)) == 29


@pytest.mark.skipif('not rcl.HAVE_PYROSETTA')
def test_db_resi_index_1coi_clash(trimerA_pose, trimerB_pose):
    # small clash_dis 2.9 necessary for this test case
    cc = PoseResnumIndex(trimerA_pose, clash_dis=2.9, nbr_atom='O')
    assert cc.contacting_resnums(trimerA_pose) is None  # clashes with self


@pytest.mark.skipif('not rcl.HAVE_PYROSETTA')
def test_db_resi_index_1coi_7(trimerA_pose, trimerB_pose):
    cc = PoseResnumIndex(trimerA_pose, clash_dis=2.9,
                         nbr_dis=7.0, nbr_atom='O')
    assert len(cc.contacting_resnums(trimerB_pose)) == 0


@pytest.mark.skipif('not rcl.HAVE_PYROSETTA')
def test_db_resi_index_1coi_765(trimerA_pose, trimerB_pose):
    cc = PoseResnumIndex(trimerA_pose, clash_dis=2.9,
                         nbr_dis=7.65, nbr_atom='O')
    assert cc.contacting_resnums(trimerB_pose) == {5, 8, 12, 15, 19}


@pytest.mark.skipif('not rcl.HAVE_PYROSETTA')
def test_db_resi_index_1coi_10(trimerA_pose, trimerB_pose):
    cc = PoseResnumIndex(trimerA_pose, clash_dis=2.9,
                         nbr_dis=10.0, nbr_atom='O')
    contacts = cc.contacting_resnums(trimerB_pose)
    assert contacts == {1, 2, 5, 8, 9, 12, 15, 16, 19, 20, 22, 23, 26}


@pytest.mark.skipif('not rcl.HAVE_PYROSETTA')
def test_db_resi_index_1coi_20(trimerA_pose, trimerB_pose):
    cc = PoseResnumIndex(trimerA_pose, clash_dis=2.9,
                         nbr_dis=20.0, nbr_atom='O')
    assert len(cc.contacting_resnums(trimerB_pose)) == 29
