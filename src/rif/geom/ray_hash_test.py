from __future__ import print_function
import rif.geom.ray_hash as rh
from pprint import pprint


def test_ray_ray_hash():
    binner4 = rh.RayToRay4dHash(1, 2, 3)
    assert str(binner4) == "RayToRay4dHash_f4i8(resl=1, lever=2, bound=3)"

    binner5 = rh.Ray5dHash(1, 2, 3)
    assert str(binner5) == "Ray5dHash_f4i8(resl=1, lever=2, bound=3)"

    binner10 = rh.RayRay10dHash(1, 2, 3)
    assert str(binner10) == "RayRay10dHash_f4i8(resl=1, lever=2, bound=3)"
    assert binner10.resl == 1
    assert binner10.lever == 2
    assert binner10.bound == 3

    for i in range(10):
        assert i == binner10.get_key(*binner10.get_center(i))  # unpack tuple
