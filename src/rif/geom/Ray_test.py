from __future__ import print_function
import numpy as np
from rif.geom import rand_ray_gaussian, Ray
import rif.eigen_types


def test_Ray_dtype():
    r = rand_ray_gaussian(np.ones(10000) * 7.0)
    n = np.linalg.norm(r['dirn']['crd'], axis=1)
    assert np.all(np.abs(1.0 - np.min(n)) < 0.0001)
    assert np.all(np.abs(1.0 - np.max(n)) < 0.0001)
    sd = np.mean(np.std(r['orig']['crd'], axis=0))
    assert abs(7 - sd) < 0.5  # generous enough, hopefully


def test_Ray():
    r = Ray(orig=[1, 2, 3], dirn=[0, 3, 0])
    assert repr(r) == "Ray(orig=[1, 2, 3], dirn=[0, 1, 0])"
    print(r.orig, r)
    assert r[0] == r.orig
    assert r[1] == r.dirn
    assert r.dirn.norm == 1
