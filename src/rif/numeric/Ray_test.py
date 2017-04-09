import numpy as np
from rif.numeric import rand_ray_gaussian


def test_Ray():
    r = rand_ray_gaussian(np.ones(10000) * 7.0)
    n = np.linalg.norm(r['dirn']['crd'], axis=1)
    assert np.all(np.abs(1.0 - np.min(n)) < 0.0001)
    assert np.all(np.abs(1.0 - np.max(n)) < 0.0001)
    sd = np.mean(np.std(r['orig']['crd'], axis=0))
    print(sd)
    assert abs(7 - sd) < 0.5  # generous enough, hopefully
