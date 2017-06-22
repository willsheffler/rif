import rif
from rif_cpp.geom import Ray
import numpy as np


def _Ray_init(self, orig=[0, 0, 0], dirn=[1, 0, 0]):
    self.orig = rif.V3(orig)
    self.dirn = rif.V3(dirn)
    self.reset_meta()


Ray.__init__ = _Ray_init


def rayorig(rays):
    assert rays['_m42']['raw'].shape[-2:] == (4, 2)
    return rays['_m42']['raw'][..., :3, 0]


def raydirn(rays):
    assert rays['_m42']['raw'].shape[-2:] == (4, 2)
    return rays['_m42']['raw'][..., :3, 1]


def ray_smalldiff(r1, r2, lever=None):
    if isinstance(r1, Ray):
        r1 = np.array([r1], dtype=Ray)
    if isinstance(r2, Ray):
        r2 = np.array([r2], dtype=Ray)
    cart = np.sqrt(np.sum((rayorig(r1) - rayorig(r2))**2, axis=1))
    # print(cart.shape)
    dots = np.einsum('ij,ij->i', raydirn(r1), raydirn(r2))
    # print(dots)
    rads = np.arccos(np.minimum(1, np.maximum(0, dots)))
    if lever is None:
        return cart, rads
    else:
        return np.sqrt(cart**2 + (rads * lever)**2)
