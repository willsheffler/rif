import sys
import pytest
import rif.dtypes
import numpy as np


def test_np_info():
    x = np.ones(5, dtype=[('crd', '(2,3)f4')])
    rif.dtypes.print_numpy_info(x)


def test_dtype_hash():
    # python2 breaks on hashing
    if sys.version_info.major is 3:
        rif.dtypes.v3f_t.__hash__()


def test_with_rif_ops():
    x = np.ones(3, dtype=rif.dtypes.v3f_t)
    with pytest.raises(TypeError):
        x + x
    with rif.dtypes.rif_ops():
        assert np.all((x + x)['crd'] == x['crd'] + x['crd'])
    with pytest.raises(TypeError):
        x * x
    assert np.all(np.arange(3) + np.arange(3) == np.arange(0, 6, 2))


def test_rif_ops_overhead():
    pass
