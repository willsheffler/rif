import rif
import rif.dtypes
import rif.geom.ray_hash as rh
from rif.geom import Ray
from pprint import pprint
import numpy as np
import pytest


def test_rayhash_invperm_error():
    h = rh.Ray5dHash(2, 10, 2000)
    # print("hsize: ", len(h))
    r = Ray(orig=[8.33408, -23.2993, -3.64781],
            dirn=[0.69189, -0.0828607, 0.717232])
    # print(r.dirn)
    # print('get key')
    i = h.get_key(r)
    # print('get center')
    c = h.get_center(i)
    # print(r.dirn)
    # print(c.dirn)
    assert c.dirn[0] > 0.6
    assert c.dirn[1] < 0.0
    assert c.dirn[2] > 0.7
