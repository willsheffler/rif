from _rif._hash._xform_hash import *
from rif import rcl
import numpy as np


class RosettaStubHash(XformHash_bt24_BCC6_X3f):
    def __init__(self, cart_resl, ang_resl, cart_bound=512):
        super().__init__(cart_resl, ang_resl, cart_bound)

    def get_key(self, x):
        return super().get_key(rcl.to_rif_stub(x))[0]

    def get_center(self, k):
        karray = np.array([k], dtype='u8')
        c = super().get_center(karray)
        return rcl.to_rosetta_stub(c)

    def __repr__(self):
        supername = 'XformHash_bt24_BCC6_X3f'
        return super().__repr__().replace(supername, 'RosettaStubHash')


class RosettaStubTorsionHash(XformAngHash_bt24_BCC6_X3f):
    def __init__(self, phi_resl, cart_resl, ang_resl, cart_bound=512):
        super().__init__(phi_resl, cart_resl, ang_resl, cart_bound)

    def get_key(self, x, ang):
        return super().get_key(rcl.to_rif_stub(x), ang)[0]

    def get_center(self, k):
        karray = np.array([k], dtype='u8')
        c = super().get_center(karray)
        return rcl.to_rosetta_stub(c[0]), c[1]

    def __repr__(self):
        supername = 'XformAngHash_bt24_BCC6_X3f'
        return super().__repr__().replace(supername, 'RosettaStubTorsionHash')


class RosettaStubTwoTorsionHash(Xform2AngHash_bt24_BCC6_X3f):
    def __init__(self, phi_resl, cart_resl, ang_resl, cart_bound=512):
        super().__init__(phi_resl, cart_resl, ang_resl, cart_bound)

    def get_key(self, x, ang1, ang2):
        return super().get_key(rcl.to_rif_stub(x), ang1, ang2)[0]

    def get_center(self, k):
        karray = np.array([k], dtype='u8')
        c = super().get_center(karray)
        return rcl.to_rosetta_stub(c[0]), c[1], c[2]

    def __repr__(self):
        supername = 'Xform2AngHash_bt24_BCC6_X3f'
        return super().__repr__().replace(supername, 'RosettaStubTwoTorsionHash')
