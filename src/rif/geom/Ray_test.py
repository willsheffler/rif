from __future__ import print_function
import numpy as np
import pytest
import rif.dtypes
from rif.geom import rand_ray_gaussian, Ray, ray_smalldiff, rayorig, raydirn
import rif.eigen_types
from numpy.testing import assert_allclose


# def test_ray_pyarray():
# a = np.array([Ray()], dtype=Ray)
# assert rif.geom.pyarray_Ray_test(a, a) == 8
# assert rif.geom.ray_hash.pyarray_Ray_test(a, a) == 7


def test_Ray_dtype():
    r = rand_ray_gaussian(sd=np.ones(10000) * 7.0)
    n = np.linalg.norm(raydirn(r), axis=1)
    assert np.all(np.abs(1.0 - np.min(n)) < 0.0001)
    assert np.all(np.abs(1.0 - np.max(n)) < 0.0001)
    sd = np.mean(np.std(rayorig(r), axis=0))
    assert abs(7 - sd) < 0.5  # generous enough, hopefully


def test_Ray():
    r = Ray(orig=[1, 2, 3], dirn=[0, 3, 0])
    # r.print_raw()
    # print(r)
    assert repr(r) == "Ray(orig=[1, 2, 3], dirn=[0, 1, 0])"
    # print(r.orig, r)
    assert r[0] == r.orig
    assert r[1] == r.dirn
    assert r.dirn.norm == 1


@pytest.mark.skipif('sys.version_info.major is 2')
def test_Ray_npy_input():
    r = Ray(orig=[1, 2, 3], dirn=[0, 3, 0])
    a = np.array([r], dtype=Ray)
    assert (rayorig(a) == [r.orig[i] for i in (0, 1, 2)]).all()
    assert (raydirn(a) == [r.dirn[i] for i in (0, 1, 2)]).all()
    assert np.all(a['_m42']['raw'][..., 3, 0] == 1)
    assert np.all(a['_m42']['raw'][..., 3, 1] == 0)


@pytest.mark.skipif('sys.version_info.major is 2')
def test_rand_ray_gaussian():
    n = 5
    sd = 10.0
    rr = rand_ray_gaussian(sd=10.0, size=n)
    assert rr.shape == (n,)
    assert_allclose(1, np.sqrt(
        np.sum(raydirn(rr)**2, axis=1)), rtol=0.00001)
    raw = rayorig(rr)
    if n > 500:
        assert all(sd - 1 < np.std(raw, axis=0))
        assert all(np.std(raw, axis=0) < sd + 1)
    cart, ang = ray_smalldiff(rr, rr)
    assert_allclose(0, cart)
    # assert_allclose(0, ang)
