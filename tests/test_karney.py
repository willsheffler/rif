from riflib.test.test_numpy import np_array_info
from riflib.sampling.orientations import read_karney_orientation_file
import riflib
import numpy as np
from riflib.sampling import orientations as ori
import pytest

def test_numpy_binding():
	np_array_info()

def test_read_karney_orientation_file():
    quat, weight = read_karney_orientation_file(
        "data/orientations/karney/c48u1.grid.gz")
    norms = np.linalg.norm(quat, axis=1)
    assert np.all(np.abs(1.0-norms) < 0.0001)
    assert np.min(weight) == np.max(weight) == 1.0

def test_karney_oddball():
    with pytest.raises(AssertionError):
        ori.quaternion_set_by_name('foobar')

def test_karney_10():
    quat, weight = riflib.sampling.orientations.quaternion_set_with_covering_radius_degrees(10.0)
    assert quat.shape == (7416,4)
    print weight[:10]
    assert abs(1.0-weight.mean()) < 0.00001

def test_karney_lengths_weights():
    # only check smallest 20 for speed
    for name in ori.karney_index.name[:20]:
        quat, weight = ori.quaternion_set_by_name(name)
        assert len(weight)
        assert abs(1.0-weight.mean()) < 1e-6
        norms = np.linalg.norm(quat, axis=1)
        assert np.all(np.abs(1.0-norms) < 0.0001)

def test_karney_by_covrad():
    assert ori.karney_name_by_radius(100) == 'c48u1'
    assert ori.karney_name_by_radius(30) == 'c48u27'
    assert ori.karney_name_by_radius(20) == 'c48u83'
    assert ori.karney_name_by_radius(10) == 'c48n309'
    assert ori.karney_name_by_radius(8) == 'c48u815'
    assert ori.karney_name_by_radius(5) == 'c48u2947'
    assert ori.karney_name_by_radius(2) == 'c48u40003'
    assert ori.karney_name_by_radius(1) == 'c48u295333'
    assert ori.karney_name_by_radius(0) == 'c48u312831'
