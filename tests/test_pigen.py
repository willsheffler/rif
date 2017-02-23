

def test_import_pigen():
    try:
        import rif.numeric.pigen
        success = True
    except:
        success = False
    assert success

from rif.numeric.pigen import Vector3f as V
from rif.numeric.pigen import Matrix33f as M
from rif.numeric.pigen import Transform3f as X


def test_Vector3f():
    u = V(1, 2, 3)
    v = V(1, 2, 3)
    assert u + v == V(2, 4, 6)
    assert u == v
    assert (u - v).norm() == 0
    u += v
    assert u == V(2, 4, 6)
    assert (u * 3).isApprox(V(6, 12, 18))


def test_Matrix33f():
    m = M()
    print(m)
    assert m == m
    assert m.isApprox(m.inverse())
    assert (m * m).isApprox(m)
    v = V(1, 2, 3)
    assert (v + v).isApprox((m + m) * v)


def test_Transform3f():
    x = X()
    print(x)
    print(x.inverse())
    assert x.isApprox(x.inverse())
    # assert 0
