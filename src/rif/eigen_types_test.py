import numpy as np
from rif.eigen_types import v3f_t, add_v3f_test


def test_eigen_types():
    a = np.zeros(10, dtype=v3f_t)
    a['x'] = np.random.randn(10)
    b = np.ones(10, dtype=v3f_t)
    assert a.shape == (10,)
    assert a['x'].shape == (10,)
    assert a['y'].shape == (10,)
    assert a['z'].shape == (10,)
    c = add_v3f_test(a, b)
    np.testing.assert_almost_equal(a['x'] + b['x'], c['x'])
    np.testing.assert_almost_equal(a['y'] + b['y'], c['y'])
    np.testing.assert_almost_equal(a['z'] + b['z'], c['z'])
