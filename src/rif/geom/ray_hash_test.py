from __future__ import print_function
import rif
import rif.dtypes
import rif.geom.ray_hash as rh
from rif.geom import Ray
from pprint import pprint
import numpy as np
import pytest


@pytest.mark.skipif('sys.version_info.major is 2')
def test_RayToRay4dHash(resl=2, n=100, lever=10):
    h = rh.RayToRay4dHash(resl, lever, bound=1000)
    for i in range(n):
        r1, r2 = rif.geom.rand_ray_gaussian(size=2, sd=10)
        r = rh.align_ray_pair(r1, r2)
        c = h.get_center(h.get_key_aligned(r))
        d = rif.geom.ray_smalldiff(r, c, lever=lever)[0]
        if d < resl:
            continue
        print("diff", d)
        print(r.dirn)
        print(c.dirn)
        print(r)
        assert d <= resl


@pytest.mark.skipif('sys.version_info.major is 2')
def test_Ray5dHash(resl=2, n=100, lever=10):
    h = rh.Ray5dHash(resl, lever, bound=1000)
    for i in range(n):
        r = rif.geom.rand_ray_gaussian(10)
        c = h.get_center(h.get_key(r))
        d = rif.geom.ray_smalldiff(r, c, lever=lever)[0]
        if d < resl:
            continue
        print("diff", d)
        print(r.dirn)
        print(c.dirn)
        print(r)
        assert d <= resl


def test_RayRay10dHash():
    h = rh.RayRay10dHash(1, 2, 3)
    assert str(h) == "RayRay10dHash_f4i8(resl=1, lever=2, bound=3)"
    assert h.resl == 1
    assert h.lever == 2
    assert h.bound == 3

# this would be better as three functions with common parts factored out???
def coverage_test_helper(cls, twoargs, N, resl, lever, sd, bound, fudge=1.0):
        # for resl in [2**i * 0.125 for i in range(6)]:
            # for lever in [2**i * 0.125 for i in range(9)]:
    h = cls(resl=resl, lever=lever, bound=bound)
    # print(cls.__name__, h.size())
    rr = rif.geom.rand_ray_gaussian(size=N, sd=sd)
    if twoargs:
        rr2 = rif.geom.rand_ray_gaussian(size=N, sd=sd)
        keys = h.get_keys(rays1=rr, rays2=rr2)
    else:
        keys = h.get_keys(rays=rr)
    cen = h.get_centers(keys)
    if twoargs and len(cen) is 2:
        diff = rif.geom.ray_smalldiff(rr, cen[0], lever=lever)
        diff2 = rif.geom.ray_smalldiff(rr2, cen[1], lever=lever)
        cen = np.where(diff > diff2, cen[0], cen[1])
        rr = np.where(diff > diff2, rr, rr2)
        diff = np.where(diff > diff2, diff, diff2)
    else:
        if twoargs:
            rr = rh.align_ray_pairs(rr, rr2)
        diff = rif.geom.ray_smalldiff(rr, cen, lever=lever)
    if np.max(diff) > resl * fudge:
        cart, rads = rif.geom.ray_smalldiff(rr, cen)
        # print('cart', cart)
        # print('rads', rads)
        # print('diff', diff)
        i = np.argmax(diff)
        print("FAIL!", cls.__name__)
        print('resl: ', resl, 'maxdiff', np.max(diff), ', cart: ',
              cart[i], ', rads: ', rads[i] * lever)
        print('ran orig', rr['orig'][i])
        print('cen orig', cen['orig'][i])
        print('ran dirn', rr['dirn'][i])
        print('cen dirn', cen['dirn'][i])
        assert np.max(diff) <= resl * fudge
    return np.max(diff)

# todo: investigate covering radius more thoroughly
def test_RayToRay4dHash_coverage():
    for resl in [2**i * .25 for i in range(4)]:
        for lever in [2**i * 0.25 for i in range(7)]:
            maxdiff = coverage_test_helper(
                rh.RayToRay4dHash, twoargs=True, resl=resl, lever=lever, sd=10,
                bound=200, N=1000, fudge=1.1)
            print('RayToRay4dHash ', resl, lever, maxdiff / resl)


def test_Ray5dHash_coverage():
    for resl in [2**i * 0.125 for i in range(6)]:
        for lever in [2**i * 0.125 for i in range(9)]:
            maxdiff = coverage_test_helper(
                rh.Ray5dHash, twoargs=False, resl=resl, lever=lever, sd=10,
                bound=200, N=1000, fudge=1.0)
            # todo: use pandas for this display!
            print('Ray5dHash ', resl, lever, maxdiff / resl)


def test_RayRay10dHash_coverage():
    for resl in [2**i * .25 for i in range(4)]:
        for lever in [2**i * 0.25 for i in range(7)]:
            maxdiff = coverage_test_helper(
                rh.RayRay10dHash, twoargs=True,
                resl=resl, lever=lever, sd=resl, bound=resl * 10,
                N=1000, fudge=1.3)
            print('RayRay10dHash ', resl, lever, maxdiff / resl)
