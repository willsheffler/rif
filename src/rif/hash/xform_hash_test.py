import pytest
from rif import rcl
from rif.hash import *
import rif.numeric.lattice
if rcl.HAVE_PYROSETTA:
    from rif.rcl import Stub


def test_XformHash_bt24_BCC6_X3f():
    xh = XformHash_bt24_BCC6_X3f(1.0, 15.0)
    assert abs(xh.cart_resl - 1.1547) < 0.001
    assert abs(xh.ori_resl - 15.0) < 0.001
    assert abs(xh.phi_resl - -1) < 0.001
    assert abs(xh.cart_bound - 512.0) < 0.001
    assert abs(xh.ori_nside - 6) < 0.001
    for i in range(100):
        x0 = xh.get_center(i)  # arbitrary xform
        k = xh.get_key(x0)
        x = xh.get_center(k)
        ktest = xh.get_key(x)
        assert ktest == k


def test_XformAngHash_bt24_BCC6_X3f():
    xh = XformAngHash_bt24_BCC6_X3f(13.0, 0.5, 10.0)
    assert abs(xh.cart_resl - 1.1547 / 2.0) < 0.001
    assert abs(xh.ori_resl - 10.0) < 0.001
    assert abs(xh.phi_resl - 13.0) < 0.001
    assert abs(xh.cart_bound - 512.0) < 0.001
    assert abs(xh.ori_nside - 9) < 0.001
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
    xh = RosettaStubHash(cart_resl=1.0, ori_resl=15.0)
    assert abs(xh.cart_resl - 1.1547) < 0.001
    assert abs(xh.ori_resl - 15.0) < 0.001
    assert abs(xh.phi_resl - -1) < 0.001
    assert abs(xh.cart_bound - 512.0) < 0.001
    assert abs(xh.ori_nside - 6) < 0.001
    assert repr(xh.grid).startswith('BCC6f4u8')
    for i in range(100):
        x0 = xh.get_center(i)  # arbitrary xform
        k = xh.get_key(x0)
        x = xh.get_center(k)
        ktest = xh.get_key(x)
        assert ktest == k


@pytest.mark.skipif('not rcl.HAVE_PYROSETTA')
def test_RosettaStubHash_with_pose(bigpose):
    xh = RosettaStubHash(cart_resl=1.0, ori_resl=15.0)
    assert abs(xh.cart_resl - 1.1547) < 0.001
    assert abs(xh.ori_resl - 15.0) < 0.001
    assert abs(xh.phi_resl - -1) < 0.001
    assert abs(xh.cart_bound - 512.0) < 0.001
    assert abs(xh.ori_nside - 6) < 0.001
    assert isinstance(xh.grid, rif.numeric.lattice.BCC6f4u8)
    assert repr(xh.grid).startswith('BCC6f4u8')
    for r in bigpose:
        stub = Stub(r.xyz('N'), r.xyz('CA'), r.xyz('C'))
        key = xh.get_key(stub)
        print(stub, key)
        # sanity check, cartesian position of bin cen should be within 1.0
        bincen = xh.get_center(key)
        assert (stub.v - bincen.v).length() < 1.0


@pytest.mark.skipif('not rcl.HAVE_PYROSETTA')
def test_RosettaStubTorsionHash():
    xh = RosettaStubTorsionHash(phi_resl=15.0, cart_resl=1.0, ori_resl=10.0)
    assert abs(xh.cart_resl - 1.1547) < 0.001
    assert abs(xh.ori_resl - 10.0) < 0.001
    assert abs(xh.phi_resl - 15.0) < 0.001
    assert abs(xh.cart_bound - 32.0) < 0.001
    assert abs(xh.ori_nside - 9) < 0.001
    assert isinstance(xh.grid, rif.numeric.lattice.BCC7f4u8)
    assert repr(xh.grid).startswith(
        'BCC7f4u8(lb=[-32 -32 -32 -0.111111 -0.111111 -0.111111 -210 ], ub=[32 32 32 1 1 1 180 ], width=[1.16364 1.16364 1.16364 0.111111 0.111111 0.111111 27.8571 ], nside=[55 55 55 10 10 10 14 ])')
    for i in range(100):
        x0, a0 = xh.get_center(i)  # arbitrary xform
        k = xh.get_key(x0, a0)
        x, a = xh.get_center(k)
        ktest = xh.get_key(x, a)
        assert ktest == k


@pytest.mark.skipif('not rcl.HAVE_PYROSETTA')
def test_RosettaStubTwoTorsionHash():
    xh = RosettaStubTwoTorsionHash(phi_resl=15.0, cart_resl=1.0, ori_resl=10.0)
    assert abs(xh.cart_resl - 1.1547) < 0.001
    assert abs(xh.ori_resl - 10.0) < 0.001
    assert abs(xh.phi_resl - 15.0) < 0.001
    assert abs(xh.cart_bound - 32.0) < 0.001
    assert abs(xh.ori_nside - 9) < 0.001
    assert isinstance(xh.grid, rif.numeric.lattice.BCC8f4u8)
    assert repr(xh.grid).startswith('BCC8f4u8(lb=[-32 -32 -32')
    for i in range(100):
        x0, a0, b0 = xh.get_center(i)  # arbitrary xform
        k = xh.get_key(x0, a0, b0)
        x, a, b = xh.get_center(k)
        ktest = xh.get_key(x, a, b)
        assert ktest == k
