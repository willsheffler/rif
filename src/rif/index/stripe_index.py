from rif import V3, Atom, _rif
import numpy as np
import functools
from rif import index


def _wrapinit(init):
    @functools.wraps(init)
    def newinit(self, radius, points, payload):
        self._init_payload = payload
        init(self, radius, points, np.arange(len(payload)))
    return newinit


def _wrappayload(func):
    def getpayload(load, x):
        if isinstance(x, tuple):
            return x[0], load[x[1]]
        else:
            return load[x]
    @functools.wraps(func)
    def newfunc(self, *args, **kwargs):
        ret = func(self, *args, **kwargs)
        return [getpayload(self._init_payload, r) for r in ret]
    return newinit


for cls in vars(_rif._index._stripe_index):
    if cls.endswith('_object'):
        print('objectifying', cls)
        cls = vars(_rif._index._stripe_index)[cls]
        cls.__init__ = _wrapinit(cls.__init__)

_typemap = {
    (V3.dtype, None): index.stripe_index_3d_V3,
    (Atom.dtype, None): index.stripe_index_3d_Atom,
    (Atom.dtype, Atom.dtype): index.stripe_index_3d_Atom,
    (V3.dtype, np.dtype('O')): index.stripe_index_3d_V3_object,
    (Atom.dtype, np.dtype('O')): index.stripe_index_3d_Atom_object,
}


def stripe_index_3d(radius, points, payload=None):
    dtype1 = np.dtype('O')
    if hasattr(points, 'dtype'):
        dtype1 = points.dtype
    dtype2 = None
    if payload is not None:
        assert len(payload) == len(points)
        dtype2 = np.dtype('O')
        if hasattr(payload, 'dtype'):
            dtype2 = payload.dtype
    return _typemap[dtype1, dtype2](radius, points, payload)
