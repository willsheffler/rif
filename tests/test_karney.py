from riflib.test.test_numpy import np_array_info
from riflib.sampling.orientations import read_karney_orientation_file
import numpy as np

def test_numpy_binding():
	np_array_info()

def test_read_karney_orientation_file():
    quats, cover = read_karney_orientation_file("data/orientations/karney/c48u1.grid.gz")
    norms = np.linalg.norm(quats,axis=1)
    assert np.all(np.abs(1.0-norms) < 0.0001)
    assert np.min(cover) == np.max(cover) == 1.0
