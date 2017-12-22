from rif.homo import *
import numpy as np
from numpy.testing import assert_allclose
import pytest


def test_homo_rotation_single():
    axis0 = np.random.randn(3)
    axis0 = axis0 / np.linalg.norm(axis0, axis=-1)
    ang0 = np.pi / 4.0
    r = hrot(list(axis0), float(ang0))
    a = fast_axis_of(r)
    n = np.linalg.norm(a, axis=-1)
    assert np.all(abs(a / n - axis0) < 0.001)
    assert np.all(abs(np.arcsin(n / 2) - ang0) < 0.001)


def test_homo_rotation_center():
    AAC = assert_allclose
    AAC([0, 2, 0, 1], hrot([1, 0, 0], 180, [0, 1, 0]) @ (0, 0, 0, 1), atol=1e-5)
    AAC([0, 1, -1, 1], hrot([1, 0, 0], 90, [0, 1, 0]) @ (0, 0, 0, 1), atol=1e-5)
    AAC([-1, 1, 2, 1], hrot([1, 1, 0], 180, [0, 1, 1]) @ (0, 0, 0, 1), atol=1e-5)
    # print(center_of(hrot([1, 0, 0], 90, [0, 1, 0])))
    # assert 0


def test_homo_rotation_array():
    shape = (1, 2, 1, 3, 4, 1, 1)
    axis0 = np.random.randn(*(shape + (3,)))
    axis0 = axis0 / np.linalg.norm(axis0, axis=-1)[..., np.newaxis]
    # todo: proper test taking neg angles into account
    #       must consider reversed axis direction
    #       also not considering close to 0 and pi/2
    ang0 = np.random.rand(*shape) * (0.99 * np.pi / 2 + 0.005 * np.pi / 2)
    r = hrot(axis0, ang0)
    # print('r', r.shape)
    a = fast_axis_of(r)
    # print('a', a.shape)
    n = np.linalg.norm(a, axis=-1)[..., np.newaxis]
    assert np.all(abs(a / n - axis0) < 0.001)
    # print('asin', np.arcsin(n / 2).shape)
    # print('ang0', ang0)
    assert np.all(abs(np.arcsin(n[..., 0] / 2) - ang0) < 0.001)


def test_htrans():
    assert htrans([1, 3, 7]).shape == (4, 4)
    assert_allclose(htrans([1, 3, 7])[:3, 3], (1, 3, 7))

    with pytest.raises(ValueError):
        htrans([4, 3, 2, 1])

    s = (2,)
    t = np.random.randn(*s, 3)
    ht = htrans(t)
    assert ht.shape == s + (4, 4)
    assert_allclose(ht[..., :3, 3], t)


def test_axis_angle_of():
    ax, an = axis_angle_of(hrot([10, 10, 0], np.pi))
    assert 1e-5 > abs(ax[0] - ax[1])
    assert 1e-5 > abs(ax[2])
    ax, an = axis_angle_of(hrot([0, 1, 0], np.pi))
    assert 1e-5 > abs(ax[0])
    assert 1e-5 > abs(ax[1]) - 1
    assert 1e-5 > abs(ax[2])

    ax, an = axis_angle_of(hrot([0, 1, 0], np.pi * 0.25))
    print(ax, an)
    assert_allclose(ax, [0, 1, 0], atol=1e-5)
    assert 1e-5 > abs(an - np.pi * 0.25)
    ax, an = axis_angle_of(hrot([0, 1, 0], np.pi * 0.75))
    print(ax, an)
    assert_allclose(ax, [0, 1, 0], atol=1e-5)
    assert 1e-5 > abs(an - np.pi * 0.75)

    ax, an = axis_angle_of(hrot([1, 0, 0], np.pi / 2))
    print(np.pi / an)
    assert 1e-5 > abs(an - np.pi / 2)


def test_axis_angle_of_rand():
    shape = (4, 5, 6, 7, 8,)
    axis = np.random.randn(*shape, 3)
    axis = axis / np.linalg.norm(axis, axis=-1)[..., np.newaxis]
    angl = np.random.random(shape) * np.pi / 2

    rot = hrot(axis, angl, dtype='f8')
    ax, an = axis_angle_of(rot)

    assert_allclose(axis, ax, rtol=1e-5)
    assert_allclose(angl, an, rtol=1e-5)
