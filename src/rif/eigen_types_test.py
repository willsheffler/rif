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
    assert np.all(a['crd'][:, 0] == 7)
    assert np.all(a['crd'][:, 1] == 8)
    assert np.all(a['crd'][:, 2] == 9)
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
    assert V3(a[0]['crd']) == b
    assert a['crd'][0, 0] == 3
    assert a['crd'][0, 1] == 4
    assert a['crd'][0, 2] == 5



@pytest.mark.skipif('sys.version_info.major is 2')
def test_eigen_v3f_dtype():
    with rif_ops():
        a = np.ones(10, dtype=v3f_t)
        a['crd'] = np.random.rand(10, 3)
        b = np.ones(10, dtype=v3f_t)
        assert a.shape == (10, )
        assert a['crd'].shape == (10, 3)
        assert_almost_equal(a['crd'] + b['crd'], (a + b)['crd'], 5)
        assert_almost_equal(a['crd'] - b['crd'], (a - b)['crd'], 5)
        assert all(np.arange(10) + np.arange(10) == np.arange(0, 20, 2))
        c = a + 2 * b
        assert_almost_equal(abs(a + 2 * b), abs(c))
        d = a[:, np.newaxis] * b
        assert np.all((2 * b)['crd'] == [2.0, 2, 2])


@pytest.mark.skipif('sys.version_info.major is 2')
def test_eigen_m3f_dtype():
    with rif_ops():
        a = np.ones(10, dtype=m3f_t)
        a['crd'] = np.random.rand(10, 3, 3)
        b = np.ones(10, dtype=m3f_t)
        assert a.shape == (10, )
        assert a['crd'].shape == (10, 3, 3)
        assert_almost_equal(a['crd'] + 2 * b['crd'], (a + 2 * b)['crd'], 5)
        assert_almost_equal(a['crd'] - 2 * b['crd'], (a - 2 * b)['crd'], 5)
        assert all(np.arange(10) + np.arange(10) == np.arange(0, 20, 2))
        c = a + b
        assert_almost_equal(abs(a + b), abs(c))
        d = a[:, np.newaxis] * b


@pytest.mark.skipif('sys.version_info.major is 2')
def test_eigen_m3f_v3f_mult():
    with rif_ops():
        a = np.empty(2, dtype=m3f_t)
        b = np.empty(2, dtype=m3f_t)

        a['crd'] = np.arange(00, 18).reshape(2, 3, 3)
        b['crd'] = np.arange(18, 36).reshape(2, 3, 3)
        assert_almost_equal(a['crd'] + b['crd'], (a + b)['crd'], 5)
        assert_almost_equal(a['crd'] - b['crd'], (a - b)['crd'], 5)
        for i in range(2):
            np_mult = a[i]['crd'].dot(b[i]['crd'])
            cp_mult = (a * b)[i]['crd']
            assert_almost_equal(np_mult, cp_mult)

        a['crd'] = np.random.randn(2, 3, 3)
        b['crd'] = np.random.randn(2, 3, 3)
        assert_almost_equal(a['crd'] + b['crd'], (a + b)['crd'], 5)
        assert_almost_equal(a['crd'] - b['crd'], (a - b)['crd'], 5)
        for i in range(2):
            np_mult = a[i]['crd'].dot(b[i]['crd'])
            cp_mult = (a * b)[i]['crd']
            assert_almost_equal(np_mult, cp_mult, 5)

        m = np.empty(2, dtype=m3f_t)
        m['crd'] = np.eye(3)
        m[1]['crd'] = m[1]['crd'] * 3
        v = np.ones(2, dtype=v3f_t)
        v2 = m * v
        assert np.all(v2[0]['crd'] == v[0]['crd'])
        assert np.all(v2[1]['crd'] == v[1]['crd'] * 3)

        v['crd'] = np.random.randn(2, 3)
        for i in range(2):
            np_mult = a[i]['crd'].dot(v[i]['crd'])
            cp_mult = (a * v)[i]['crd']
            assert_almost_equal(np_mult, cp_mult, 5)


def test_eigen_fields():
    x = np.zeros(4, dtype=rif.eigen_types.test_t)
    assert x.shape == (4,)
    assert x['a']['crd'].shape == (4, 3)
    assert x['b']['crd'].shape == (4, 3)
    assert x['i'].dtype == 'i4'
    assert x['f'].dtype == 'f4'
