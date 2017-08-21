from rif.legacy import xyzMath
from rif.legacy.xyzMath import Vec
import doctest


def test_xyzMath():
    v = Vec(1, 0, 0)
    assert v + 2 * v == Vec(3, 0, 0)


def test_xyzMath_doctest():
    result = doctest.testmod(xyzMath)
    assert result.failed == 0
