import numpy as np
from rif.test.dtype_test import struct_t


def test_eigen_types():
    x = np.zeros(10, dtype=struct_t)
    assert x.shape == (10,)
    assert x['FIELD1'].shape == (10,)
    assert x['FIELD2'].shape == (10,)
