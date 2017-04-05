import numpy as np
from rif.eigen_types import fxyz_t


def test_eigen_types():
    x = np.zeros(10, dtype=fxyz_t)
    assert x.shape == (10,)
    assert x['x'].shape == (10,)
