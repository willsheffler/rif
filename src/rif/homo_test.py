from rif.homo import *
import numpy as np
from numpy.testing import assert_array_almost_equal as assert_aae


def test_homo_rotation_single():
    axis0 = np.random.randn(3)
    axis0 = axis0 / np.linalg.norm(axis0, axis=-1)
    ang0 = np.pi / 2.0
    r = homo_rotation(list(axis0), float(ang0))
    a = axis_of(r)
    n = np.linalg.norm(a, axis=-1)
    assert np.all(abs(a / n - axis0) < 0.001)
    assert np.all(abs(np.arcsin(n / 2) - ang0) < 0.001)


def test_homo_rotation_center():
    HR = homo_rotation
    assert_aae([0, 2, 0, 1], HR([1, 0, 0], 180, [0, 1, 0]) @ (0, 0, 0, 1))
    assert_aae([0, 1, -1, 1], HR([1, 0, 0], 90, [0, 1, 0]) @ (0, 0, 0, 1))
    assert_aae([-1, 1, 2, 1], HR([1, 1, 0], 180, [0, 1, 1]) @ (0, 0, 0, 1))


def test_homo_rotation_array():
    shape = (1, 2, 1, 3, 4, 1, 1)
    axis0 = np.random.randn(*(shape + (3,)))
    axis0 = axis0 / np.linalg.norm(axis0, axis=-1)[..., np.newaxis]
    # todo: proper test taking neg angles into account
    #       must consider reversed axis direction
    #       also not considering close to 0 and pi/2
    ang0 = np.random.rand(*shape) * (0.99 * np.pi / 2 + 0.005 * np.pi / 2)
    r = homo_rotation(axis0, ang0)
    # print('r', r.shape)
    a = axis_of(r)
    # print('a', a.shape)
    n = np.linalg.norm(a, axis=-1)[..., np.newaxis]
    assert np.all(abs(a / n - axis0) < 0.001)
    # print('asin', np.arcsin(n / 2).shape)
    # print('ang0', ang0)
    assert np.all(abs(np.arcsin(n[..., 0] / 2) - ang0) < 0.001)
