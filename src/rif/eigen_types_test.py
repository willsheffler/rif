from __future__ import print_function
import pytest
import numpy as np
import rif
from rif.dtypes import rif_ops
from rif.eigen_types import *
from pprint import pprint
from numpy.testing import assert_almost_equal


def test_V3_numpy():
    # print([V3()] * 5)
    a = np.array([(V3((7, 8, 9)),)] * 4, dtype=V3)
    assert np.all(a['raw'][:, 0] == 7)
    assert np.all(a['raw'][:, 1] == 8)
    assert np.all(a['raw'][:, 2] == 9)
    # print(a.dtype)
    # print(a[0])
    # assert repr(a[0]) == ''

    print(type(a[0]))
    print(type(np.asscalar(a[0])))
    assert V3(a[2][0]) == a[2]

    u = V3([1, 2, 3])
    v = V3([1, 2, 3])
    assert u == v
    assert len(v) == 3
    assert V3.dtype == v3d_t
    assert V3().dtype == v3d_t

    # assert 0


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


@pytest.mark.skipif('sys.version_info.major is 2')
def test_eigen_v3f_dtype():
    with rif_ops():
        a = np.ones(10, dtype=v3f_t)
        a['raw'] = np.random.rand(10, 3)
        b = np.ones(10, dtype=v3f_t)
        assert a.shape == (10, )
        assert a['raw'].shape == (10, 3)
        assert_almost_equal(a['raw'] + b['raw'], (a + b)['raw'], 5)
        assert_almost_equal(a['raw'] - b['raw'], (a - b)['raw'], 5)
        assert all(np.arange(10) + np.arange(10) == np.arange(0, 20, 2))
        c = a + 2 * b
        assert_almost_equal(abs(a + 2 * b), abs(c))
        d = a[:, np.newaxis] * b
        assert np.all((2 * b)['raw'] == [2.0, 2, 2])


@pytest.mark.skipif('sys.version_info.major is 2')
def test_eigen_m3f_dtype():
    with rif_ops():
        a = np.ones(10, dtype=m3f_t)
        a['raw'] = np.random.rand(10, 3, 3)
        b = np.ones(10, dtype=m3f_t)
        assert a.shape == (10, )
        assert a['raw'].shape == (10, 3, 3)
        assert_almost_equal(a['raw'] + 2 * b['raw'], (a + 2 * b)['raw'], 5)
        assert_almost_equal(a['raw'] - 2 * b['raw'], (a - 2 * b)['raw'], 5)
        assert all(np.arange(10) + np.arange(10) == np.arange(0, 20, 2))
        c = a + b
        assert_almost_equal(abs(a + b), abs(c))
        d = a[:, np.newaxis] * b


@pytest.mark.skipif('sys.version_info.major is 2')
def test_eigen_m3f_v3f_mult():
    with rif_ops():
        a = np.empty(2, dtype=m3f_t)
        b = np.empty(2, dtype=m3f_t)

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

        m = np.empty(2, dtype=m3f_t)
        m['raw'] = np.eye(3)
        m[1]['raw'] = m[1]['raw'] * 3
        v = np.ones(2, dtype=v3f_t)
        v2 = m * v
        assert np.all(v2[0]['raw'] == v[0]['raw'])
        assert np.all(v2[1]['raw'] == v[1]['raw'] * 3)

        v['raw'] = np.random.randn(2, 3)
        for i in range(2):
            np_mult = a[i]['raw'].dot(v[i]['raw'])
            cp_mult = (a * v)[i]['raw']
            assert_almost_equal(np_mult, cp_mult, 5)


def test_eigen_fields():
    x = np.zeros(4, dtype=rif.eigen_types.test_t)
    assert x.shape == (4,)
    assert x['a']['raw'].shape == (4, 3)
    assert x['b']['raw'].shape == (4, 3)
    assert x['i'].dtype == 'i4'
    assert x['f'].dtype == 'f4'
