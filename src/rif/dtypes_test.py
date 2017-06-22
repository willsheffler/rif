import sys
import pytest
import rif
from rif import V3, M3, X3
from rif.eigen_types import Xc3
from rif.dtypes import RifOperators, RifOperatorsDisabled
import numpy as np
from numpy.testing import assert_almost_equal


def test_plus_equals_forwarding():
    try:
        a = np.array([1, 2, 3, 4], dtype='int64')
        a += a
    except TypeError as e:
        pytest.fail('unexpected TypeError', e)


def test_dtype_hash():
    # python2 breaks on hashing
    if sys.version_info.major is 3:
        V3.dtype.__hash__()


def test_V3_numpy():
    # print([V3()] * 5)
    a = np.array([(V3((7, 8, 9)),)] * 4, dtype=V3)
    assert np.all(a['raw'][:, 0] == 7)
    assert np.all(a['raw'][:, 1] == 8)
    assert np.all(a['raw'][:, 2] == 9)

    if sys.version_info.major is 3:
        assert np.all(abs(a) == np.linalg.norm(a['raw'], axis=1))

    # print(a.dtype)
    # print(a[0])
    # assert repr(a[0]) == ''

    # print(type(a[0]))
    # print(type(np.asscalar(a[0])))
    assert V3(a[2][0]) == a[2]

    u = V3([1, 2, 3])
    v = V3([1, 2, 3])
    assert u == v
    assert len(v) == 3


def eigen_V3_test_helper():
    a = np.ones(10, dtype=V3)
    a['raw'] = np.random.rand(10, 3)
    b = np.ones(10, dtype=V3)
    assert a.shape == (10, )
    assert a['raw'].shape == (10, 3)
    assert_almost_equal(a['raw'] + b['raw'], (a + b)['raw'], 5)
    assert_almost_equal(a['raw'] - b['raw'], (a - b)['raw'], 5)
    assert all(np.arange(10) + np.arange(10) == np.arange(0, 20, 2))
    c = a + 2 * b
    assert_almost_equal(abs(a + 2.0 * b), abs(c))
    # d = a[:, np.newaxis] * b
    assert np.all((2.0 * b)['raw'] == [2.0, 2, 2])


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


def test_X3_numpy():
    n = 2
    a = np.zeros(n, dtype=X3)
    assert a['raw'].shape == (n, 4, 4)
    a = np.zeros(n, dtype=Xc3)
    assert a['raw'].shape == (n, 3, 4)


def test_numpy_delegation():
    assert all((3 * np.arange(3)) == np.array([0, 3, 6]))


@pytest.mark.skipif('sys.version_info.major is 2')
def test_V3_numpy_assign():
    # print([V3()] * 5)
    a = np.array([((7, 8, 9),)] * 4, dtype=V3)
    b = V3([3, 4, 5])
    a[0] = b
    assert V3(a[0]['raw']) == b
    assert a['raw'][0, 0] == 3
    assert a['raw'][0, 1] == 4
    assert a['raw'][0, 2] == 5


def test_eigen_V3_dtype():
    with rif.dtypes.RifOperators():
        eigen_V3_test_helper()


def test_eigen_M3_dtype():
    with RifOperators():
        a = np.ones(10, dtype=M3)
        a['raw'] = np.random.rand(10, 3, 3)
        b = np.ones(10, dtype=M3)
        assert a.shape == (10, )
        assert a['raw'].shape == (10, 3, 3)
        assert_almost_equal(a['raw'] + 2 * b['raw'], (a + 2 * b)['raw'], 5)
        assert_almost_equal(a['raw'] - 2 * b['raw'], (a - 2 * b)['raw'], 5)
        assert all(np.arange(10) + np.arange(10) == np.arange(0, 20, 2))
        c = a + b
        assert_almost_equal(abs(a + b), abs(c))
        # d = a[:, np.newaxis] * b


def test_eigen_M3_V3_mult():
    with RifOperators():
        a = np.empty(2, dtype=M3)
        b = np.empty(2, dtype=M3)

        a['raw'] = np.arange(00, 18).reshape(2, 3, 3)
        b['raw'] = np.arange(18, 36).reshape(2, 3, 3)
        assert_almost_equal(a['raw'] + b['raw'], (a + b)['raw'], 5)
        assert_almost_equal(a['raw'] - b['raw'], (a - b)['raw'], 5)
        for i in range(2):
            np_mult = a[i]['raw'].dot(b[i]['raw'])
            cp_mult = (a * b)[i]['raw']
            assert_almost_equal(np_mult, cp_mult)

        a['raw'] = np.random.randn(2, 3, 3)
        b['raw'] = np.random.randn(2, 3, 3)
        assert_almost_equal(a['raw'] + b['raw'], (a + b)['raw'], 5)
        assert_almost_equal(a['raw'] - b['raw'], (a - b)['raw'], 5)
        for i in range(2):
            np_mult = a[i]['raw'].dot(b[i]['raw'])
            cp_mult = (a * b)[i]['raw']
            assert_almost_equal(np_mult, cp_mult, 5)

        m = np.empty(2, dtype=M3)
        m['raw'] = np.eye(3)
        m[1]['raw'] = m[1]['raw'] * 3
        v = np.ones(2, dtype=V3)
        v2 = m * v
        assert np.all(v2[0]['raw'] == v[0]['raw'])
        assert np.all(v2[1]['raw'] == v[1]['raw'] * 3)

        v['raw'] = np.random.randn(2, 3)
        for i in range(2):
            np_mult = a[i]['raw'].dot(v[i]['raw'])
            cp_mult = (a * v)[i]['raw']
            assert_almost_equal(np_mult, cp_mult, 5)
