from rif.hash import *
from rif import rcl
import pytest


def test_XformHash_bt24_BCC6_X3f():
    xh = XformHash_bt24_BCC6_X3f(1.0, 15.0)
    name = 'XformHash_bt24_BCC6_X3f(cart_resl=1, ang_resl=15, cart_bound=512)'
    assert repr(xh) == name
    for i in range(100):
        x0 = xh.get_center(i)  # arbitrary xform
        k = xh.get_key(x0)
        x = xh.get_center(k)
        ktest = xh.get_key(x)
        assert ktest == k


@pytest.mark.skipif('not rcl.HAVE_PYROSETTA')
def test_XformHash_bt24_BCC6_Rosetta():
    xh = XformHash_bt24_BCC6_Rosetta(1.0, 15.0)
    name = 'XformHash_bt24_BCC6_Rosetta(cart_resl=1, ang_resl=15, cart_bound=512)'
    assert repr(xh) == name
    for i in range(100):
        x0 = xh.get_center(i)  # arbitrary xform
        k = xh.get_key(x0)
        x = xh.get_center(k)
        ktest = xh.get_key(x)
        assert ktest == k
