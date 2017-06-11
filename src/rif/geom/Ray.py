from rif_cpp.geom import Ray
import numpy as np


def _Ray_init(self, orig=[0, 0, 0], dirn=[1, 0, 0]):
    for i in (0, 1, 2):
        self.orig[i] = orig[i]
        self.dirn[i] = dirn[i]
    self.dirn.normalize()


Ray.__init__ = _Ray_init


def ray_smalldiff(r1, r2, lever=None):
    if isinstance(r1, Ray):
        print(type(r1.orig))
        r1 = np.array([(r1.orig, r1.dirn)], dtype=Ray)
    if isinstance(r2, Ray):
        r2 = np.array([(r2.orig, r2.dirn)], dtype=Ray)
    cart = np.sqrt(np.sum((r1['orig']['raw'] - r2['orig']['raw'])**2, axis=1))
    # print(cart.shape)
    dots = np.einsum('ij,ij->i', r1['dirn']['raw'], r2['dirn']['raw'])
    # print(dots)
    rads = np.arccos(np.minimum(1, np.maximum(0, dots)))
    if lever is None:
        return cart, rads
    else:
        return np.sqrt(cart**2 + (rads * lever)**2)
