import numpy as np
import rif
from rif.eigen_types import V3, X3, x3_identity, x3_inverse
from rif.dtypes import RifOperators


def test_np_info():
    x = np.ones(5, dtype=[('raw', '(2,3)f4')])
    rif.eigen_types.print_numpy_info(x)


def test_V3():
    assert V3(1, 2, 3) == V3([1, 2, 3])
    assert V3(1, 2, 3) != V3([1, 2, 4])


def test_inverse():
    with RifOperators():
        assert x3_inverse(x3_identity) == x3_identity
        x = np.array([1, 0, 0, 10, 0, 0, 1, 11, 0, 1, 0, 12,
                      0, 0, 0, 1] * 10, dtype='f').view(X3)
        assert (x * x3_inverse(x) == x3_identity).all()
        assert (x3_inverse(x) * x == x3_identity).all()
