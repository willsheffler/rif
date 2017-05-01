from rif.geom.rand_geom import *
import numpy as np


def test_rand_normal():
    lens = np.ones(100)
    vecs = rand_normal(lens)
    assert np.all(np.abs(np.linalg.norm(vecs['crd'], axis=1) - lens) < 0.0001)
