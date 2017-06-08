import sys
import pytest
from rif import V3, dtypes as dt
from rif.dtypes import RifOperators, RifOperatorsDisabled
import numpy as np


def test_dtype_hash():
    # python2 breaks on hashing
    if sys.version_info.major is 3:
        V3.dtype.__hash__()


def test_with_rif_ops():
    with RifOperatorsDisabled():
        x = np.ones(3, dtype=V3)
        with pytest.raises(TypeError):
            x + x
        with RifOperators():
            assert np.all((x + x)['raw'] == x['raw'] + x['raw'])
        with pytest.raises(TypeError):
            x * x
        assert np.all(np.arange(3) + np.arange(3) == np.arange(0, 6, 2))


def test_rif_ops_overhead():
    n = 2
    np.zeros(n, dtype=V3)
    x = np.random.uniform(size=100).reshape((10, 10))
    print(x)

    # assert 0
