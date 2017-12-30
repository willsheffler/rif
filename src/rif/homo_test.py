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


def test_is_valid_rays():
    assert not is_valid_rays([[0, 0, 0, 0], [1, 0, 0, 0]])
    assert not is_valid_rays([[0, 0, 0, 1], [0, 0, 0, 0]])
    assert not is_valid_rays([[0, 0, 0, 1], [0, 3, 0, 0]])
    assert is_valid_rays([[0, 0, 0, 1], [0, 1, 0, 0]])


def test_random_rays():
    r = random_rays()
    assert np.all(r[..., :, 3] == (1, 0))
    assert r.shape == (2, 4)
    assert_allclose(np.linalg.norm(r[..., 1, :3], axis=-1), 1)

    r = random_rays(shape=(5, 6, 7))
    assert np.all(r[..., :, 3] == (1, 0))
    assert r.shape == (5, 6, 7, 2, 4)
    assert_allclose(np.linalg.norm(r[..., 1, :3], axis=-1), 1)


def test_proj_prep():
    assert_allclose([2, 3, 0], proj_perp([0, 0, 1], [2, 3, 99]))
    assert_allclose([2, 3, 0], proj_perp([0, 0, 2], [2, 3, 99]))
    a, b = np.random.randn(2, 5, 6, 7, 3)
    pp = proj_perp(a, b)
    assert_allclose(hdot(a, pp), 0, atol=1e-5)


def test_point_in_plane():
    plane = random_rays((5, 6, 7))
    assert np.all(point_in_plane(plane, plane[..., 0, :3]))
    pt = proj_perp(plane[..., 1, :3], np.random.randn(3))
    assert np.all(point_in_plane(plane, plane[..., 0, :3] + pt))


def test_ray_in_plane():
    plane = random_rays((5, 6, 7))
    dirn = proj_perp(plane[..., 1, :3], np.random.randn(5, 6, 7, 3))
    ray = hray(plane[..., 0, :3] + np.cross(plane[..., 1, :3], dirn) * 7, dirn)
    assert np.all(ray_in_plane(plane, ray))


def test_intersect_planes():
    with pytest.raises(ValueError):
        intersect_planes(np.array([[0, 0, 0, 2], [0, 0, 0, 0]]),
                         np.array([[0, 0, 0, 1], [0, 0, 0, 0]]))
    with pytest.raises(ValueError):
        intersect_planes(np.array([[0, 0, 0, 1], [0, 0, 0, 0]]),
                         np.array([[0, 0, 0, 1], [0, 0, 0, 1]]))
    with pytest.raises(ValueError):
        intersect_planes(np.array([[0, 0, 1], [0, 0, 0, 0]]),
                         np.array([[0, 0, 1], [0, 0, 0, 1]]))
    with pytest.raises(ValueError):
        intersect_planes(np.array(9 * [[[0, 0, 0, 1], [0, 0, 0, 0]]]),
                         np.array(2 * [[[0, 0, 0, 1], [0, 0, 0, 0]]]))

    # isct, sts = intersect_planes(np.array(9 * [[[0, 0, 0, 1], [1, 0, 0, 0]]]),
        # np.array(9 * [[[0, 0, 0, 1], [1, 0, 0, 0]]]))
    # assert isct.shape[:-2] == sts.shape == (9,)
    # assert np.all(sts == 2)

    # isct, sts = intersect_planes(np.array([[1, 0, 0, 1], [1, 0, 0, 0]]),
        # np.array([[0, 0, 0, 1], [1, 0, 0, 0]]))
    # assert sts == 1

    isct, sts = intersect_planes(np.array([[0, 0, 0, 1], [1, 0, 0, 0]]),
                                 np.array([[0, 0, 0, 1], [0, 1, 0, 0]]))
    assert sts == 0
    assert isct[0, 2] == 0
    assert np.all(abs(isct[1, :3]) == (0, 0, 1))

    isct, sts = intersect_planes(np.array([[0, 0, 0, 1], [1, 0, 0, 0]]),
                                 np.array([[0, 0, 0, 1], [0, 0, 1, 0]]))
    assert sts == 0
    assert isct[0, 1] == 0
    assert np.all(abs(isct[1, :3]) == (0, 1, 0))

    isct, sts = intersect_planes(np.array([[0, 0, 0, 1], [0, 1, 0, 0]]),
                                 np.array([[0, 0, 0, 1], [0, 0, 1, 0]]))
    assert sts == 0
    assert isct[0, 0] == 0
    assert np.all(abs(isct[1, :3]) == (1, 0, 0))

    isct, sts = intersect_planes(np.array([[7, 0, 0, 1], [1, 0, 0, 0]]),
                                 np.array([[0, 9, 0, 1], [0, 1, 0, 0]]))
    assert sts == 0
    assert_allclose(isct[0, :3], [7, 9, 0])
    assert_allclose(abs(isct[1, :3]), [0, 0, 1])

    isct, sts = intersect_planes(np.array([[0, 0, 0, 1], hnormalized([1, 1, 0, 0])]),
                                 np.array([[0, 0, 0, 1], hnormalized([0, 1, 1, 0])]))
    assert sts == 0
    assert_allclose(abs(isct[1, :3]), hnormalized([1, 1, 1]))

    p1 = np.array([[2, 0, 0, 1], hnormalized([1, 0, 0, 0])])
    p2 = np.array([[0, 0, 0, 1], hnormalized([0, 0, 1, 0])])
    isct, sts = intersect_planes(p1, p2)
    assert sts == 0
    assert np.all(ray_in_plane(p1, isct))
    assert np.all(ray_in_plane(p2, isct))

    p1 = np.array([[0.39263901, 0.57934885, -0.7693232, 1.],
                   [-0.80966465, -0.18557869, 0.55677976, 0.]])
    p2 = np.array([[0.14790894, -1.333329, 0.45396509, 1.],
                   [-0.92436319, -0.0221499, 0.38087016, 0.]])
    isct, sts = intersect_planes(p1, p2)
    assert sts == 0
    assert np.all(ray_in_plane(p1, isct))
    assert np.all(ray_in_plane(p2, isct))


def test_intersect_planes_rand():
    plane1, plane2 = random_rays(shape=(2, 1))
    plane1[..., 0, :3] = 0
    plane2[..., 0, :3] = 0
    isect, status = intersect_planes(plane1, plane2)
    assert np.all(status == 0)
    assert np.all(ray_in_plane(plane1, isect))
    assert np.all(ray_in_plane(plane2, isect))

    plane1, plane2 = random_rays(shape=(2, 1))
    plane1[..., 1, :3] = hnormalized([0, 0, 1])
    plane2[..., 1, :3] = hnormalized([0, 1, 0])
    isect, status = intersect_planes(plane1, plane2)
    assert np.all(status == 0)
    assert np.all(ray_in_plane(plane1, isect))
    assert np.all(ray_in_plane(plane2, isect))

    plane1, plane2 = random_rays(shape=(2, 5, 6, 7, 8, 9))
    isect, status = intersect_planes(plane1, plane2)
    assert np.all(status == 0)
    assert np.all(ray_in_plane(plane1, isect))
    assert np.all(ray_in_plane(plane2, isect))
