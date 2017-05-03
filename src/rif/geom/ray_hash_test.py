from __future__ import print_function
import rif.geom.ray_hash as rh
from pprint import pprint


def test_ray_ray_hash():
    binner4 = rh.RayToRayBinner4D(1, 2, 3)
    assert str(binner4) == "RayToRayBinner4D_f8i8(resl=1, lever=2, bound=3)"

    binner5 = rh.RayBinner5D(1, 2, 3)
    assert str(binner5) == "RayBinner5D_f8i8(resl=1, lever=2, bound=3)"

    binner10 = rh.RayRayBinner10D(1, 2, 3)
    assert str(binner10) == "RayRayBinner10D_f8i8(resl=1, lever=2, bound=3)"
    assert binner10.resl == 1
    assert binner10.lever == 2
    assert binner10.bound == 3

    for i in range(10):
        assert i == binner10.get_key(binner10.get_center(i))
