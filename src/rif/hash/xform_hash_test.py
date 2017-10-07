import pytest
from rif import rcl
from rif.hash import *
if rcl.HAVE_PYROSETTA:
    from rif.rcl import Stub


def test_XformHash_bt24_BCC6_X3f():
    xh = XformHash_bt24_BCC6_X3f(1.0, 15.0)
    for i in range(100):
        x0 = xh.get_center(i)  # arbitrary xform
        k = xh.get_key(x0)
        x = xh.get_center(k)
        ktest = xh.get_key(x)
        assert ktest == k


def test_XformAngHash_bt24_BCC6_X3f():
    xh = XformAngHash_bt24_BCC6_X3f(15.0, 1.0, 15.0)
    for i in range(100):
        x0 = xh.get_center(i)  # arbitrary xform
        print(type(x0[1]))
        # k = xh.get_key(x0[0], x0[1])
        k = xh.get_key(x0)
        x = xh.get_center(k)
        # ktest = xh.get_key(x[0], x[1])
        ktest = xh.get_key(x)
        assert ktest == k


def test_Xform2AngHash_bt24_BCC6_X3f():
    xh = Xform2AngHash_bt24_BCC6_X3f(15.0, 1.0, 15.0)
    for i in range(100):
        x0 = xh.get_center(i)  # arbitrary xform
        k = xh.get_key(x0)
        x = xh.get_center(k)
        ktest = xh.get_key(x)
        assert ktest == k


@pytest.mark.skipif('not rcl.HAVE_PYROSETTA')
def test_RosettaStubHash():
    xh = RosettaStubHash(cart_resl=1.0, ang_resl=15.0)
    for i in range(100):
        x0 = xh.get_center(i)  # arbitrary xform
        k = xh.get_key(x0)
        x = xh.get_center(k)
        ktest = xh.get_key(x)
        assert ktest == k


@pytest.mark.skipif('not rcl.HAVE_PYROSETTA')
def test_RosettaStubHash_with_pose(bigpose):
    stub_hasher = RosettaStubHash(cart_resl=1.0, ang_resl=15.0)
    for r in bigpose:
        stub = Stub(r.xyz('N'), r.xyz('CA'), r.xyz('C'))
        key = stub_hasher.get_key(stub)
        print(stub, key)
        # sanity check, cartesian position of bin cen should be within 1.0
        bincen = stub_hasher.get_center(key)
        assert (stub.v - bincen.v).length() < 1.0


@pytest.mark.skipif('not rcl.HAVE_PYROSETTA')
def test_RosettaStubTorsionHash():
    xh = RosettaStubTorsionHash(phi_resl=15.0, cart_resl=1.0, ang_resl=10.0)
    for i in range(100):
        x0, a0 = xh.get_center(i)  # arbitrary xform
        k = xh.get_key(x0, a0)
        x, a = xh.get_center(k)
        ktest = xh.get_key(x, a)
        assert ktest == k


@pytest.mark.skipif('not rcl.HAVE_PYROSETTA')
def test_RosettaStubTwoTorsionHash():
    xh = RosettaStubTwoTorsionHash(phi_resl=15.0, cart_resl=1.0, ang_resl=10.0)
    for i in range(100):
        x0, a0, b0 = xh.get_center(i)  # arbitrary xform
        k = xh.get_key(x0, a0, b0)
        x, a, b = xh.get_center(k)
        ktest = xh.get_key(x, a, b)
        assert ktest == k
