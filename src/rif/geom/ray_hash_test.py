from __future__ import print_function
import rif.geom.ray_hash as rh
from pprint import pprint


def test_ray_ray_hash():
    binner4 = rh.RayToRayHash4D(1, 2, 3)
    assert str(binner4) == "RayToRayHash4D_f4i8(resl=1, lever=2, bound=3)"

    binner5 = rh.RayHash5D(1, 2, 3)
    assert str(binner5) == "RayHash5D_f4i8(resl=1, lever=2, bound=3)"

    binner10 = rh.RayRayHash10D(1, 2, 3)
    assert str(binner10) == "RayRayHash10D_f4i8(resl=1, lever=2, bound=3)"
    assert binner10.resl == 1
    assert binner10.lever == 2
    assert binner10.bound == 3

    for i in range(10):
        assert i == binner10.get_key(binner10.get_center(i))
