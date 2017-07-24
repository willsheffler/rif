import rif
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


def _wrappayload(cls, funcname):
    orig = getattr(cls, funcname)
    def dowrap(func):
        def getpayload(load, x):
            if isinstance(x, tuple):
                return x[0], load[x[1]]
            else:
                return load[x]
        @functools.wraps(func)
        def newfunc(self, *args, **kwargs):
            ret = func(self, *args, **kwargs)
            if hasattr(ret, '__iter__'):
                return [getpayload(self._init_payload, r) for r in ret]
            else:
                return getpayload(self._init_payload, ret)
        return newfunc
    setattr(cls, funcname, dowrap(orig))


_typemap_3d = dict()


def get_dtype(name):
    if name == 'object':
        return np.dtype('O')
    else:
        return vars(rif)[name].dtype

for clsname in vars(_rif._index._stripe_index):
    cls = vars(_rif._index._stripe_index)[clsname]
    if clsname.endswith('_object'):
        cls.__init__ = _wrapinit(cls.__init__)
        _wrappayload(cls, '_raw_payload')
        _wrappayload(cls, 'neighbors')
    if clsname.startswith('stripe_index_3d_'):
        dtypes = tuple([get_dtype(x) for x in clsname.split('_')[3:]])
        if len(dtypes) is 1:
            dtypes = (dtypes[0], None)
        _typemap_3d[dtypes] = cls


def convert_to_arrayV3(iterable):
    seq = list(iterable)
    array = np.zeros((len(seq), 3), dtype='f4')
    for i, v in enumerate(seq):
        array[i, 0] = v[0]
        array[i, 1] = v[1]
        array[i, 2] = v[2]
    return array.view(V3)


def stripe_index_3d(radius, points, payload=None):
    """create a 3D stripe index of points[payload] based on dtype of inputs

    stripe index is a fast index of 3D coords (possible enhanced with metadata)
    which can quickly get all the neighbors within a certain radius, which must
    be specified at constructin

    Arguments:
        radius {float} -- radius of neighbor lookup
        points {np.ndarray} -- array of dtype V3, Atom, etc...
        payload {np.ndarray} -- array of any dtype, or sequence of arbitrary
                                objects (default: {No payload})

    Returns:
        wrapped c++ stripe index of appropriate type, depending on dtype of
        input arrays

    Raises:
        ValueError -- if points or payload types not understood
    """
    dtype1 = np.dtype('O')
    if hasattr(points, 'dtype'):
        dtype1 = points.dtype
    dtype2 = None
    if payload is not None:
        assert len(payload) == len(points)
        dtype2 = np.dtype('O')
        if hasattr(payload, 'dtype'):
            dtype2 = payload.dtype
        else:
            payload = np.array(payload, dtype='O')
            if not hasattr(points, 'dtype'):
                # if payload is sequence of objects, allow points
                # to be an iterable of sequences of len at least 3
                # for example, a list of rosetta xyzVectors
                points = convert_to_arrayV3(points)
                dtype1 = V3.dtype
    try:
        return _typemap_3d[dtype1, dtype2](radius, points, payload)
    except KeyError:

        raise ValueError(
            "stripe_index_3d only supports numpy arrays (can have " +
            "fancy dtypes) or object payload (which is inefficient), not: " +
            str(type(points)) + ', ' + str(type(payload)))
