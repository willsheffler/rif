from __future__ import print_function
import numpy as np
import pytest
import rif.dtypes
from rif.geom import rand_ray_gaussian, Ray, ray_smalldiff
import rif.eigen_types
from numpy.testing import assert_allclose


def test_Ray_dtype():
    r = rand_ray_gaussian(sd=np.ones(10000) * 7.0)
    n = np.linalg.norm(r['dirn']['raw'], axis=1)
    assert np.all(np.abs(1.0 - np.min(n)) < 0.0001)
    assert np.all(np.abs(1.0 - np.max(n)) < 0.0001)
    sd = np.mean(np.std(r['orig']['raw'], axis=0))
    assert abs(7 - sd) < 0.5  # generous enough, hopefully


def test_Ray():
    r = Ray(orig=[1, 2, 3], dirn=[0, 3, 0])
    assert repr(r) == "Ray(orig=[1, 2, 3], dirn=[0, 1, 0])"
    print(r.orig, r)
    assert r[0] == r.orig
    assert r[1] == r.dirn
    assert r.dirn.norm == 1


@pytest.mark.skipif('sys.version_info.major is 2')
def test_rand_ray_gaussian():
    n = 5
    sd = 10.0
    rr = rand_ray_gaussian(sd=10.0, size=n)
    assert rr.shape == (n,)
    assert_allclose(1, np.sqrt(
        np.sum(rr['dirn']['raw']**2, axis=1)), rtol=0.00001)
    raw = rr['orig']['raw']
    if n > 500:
        assert all(sd - 1 < np.std(raw, axis=0))
        assert all(np.std(raw, axis=0) < sd + 1)
    cart, ang = ray_smalldiff(rr, rr)
    assert_allclose(0, cart)
    # assert_allclose(0, ang)
