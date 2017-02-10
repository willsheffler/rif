


def test_import_pigen():
    try:
        import pigen
        success = True
    except:
        success = False
    assert success

from pigen import Vector3f as V
from pigen import Matrix33f as M
from pigen import Transform3f as X

def test_Vector3f():
    u = V(1,2,3)
    v = V(1,2,3)
    print u
    assert u == v


def test_Matrix33f():
    m = M()
    print m
    assert m == m
    # assert 0

def test_Transform3f():
    x = X()
    print x
    print x.inverse()
    assert x.isApprox(x.inverse())
    # assert 0