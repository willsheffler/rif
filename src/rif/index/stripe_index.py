import rif
from rif import V3, Atom, _rif
import numpy as np
import functools


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
            # payload = np.array(payload, dtype='O')
            if not hasattr(points, 'dtype'):
                # if payload is sequence of objects, allow points
                # to be an iterable of 'index-ables' of len at least 3
                # for example, a list of rosetta xyzVectors
                points = _convert_to_arrayV3(points)
                dtype1 = V3.dtype
    try:
        return _typemap_3d[dtype1, dtype2](radius, points, payload)
    except KeyError:
        raise ValueError(
            "stripe_index_3d only supports numpy arrays (can have " +
            "fancy dtypes) or object payload (which is inefficient), not: " +
            str(type(points)) + ', ' + str(type(payload)))


def _enable_init_pyobject_store(init):
    @functools.wraps(init)
    def newinit(self, radius, points, payload):
        self._init_payload = list(payload)
        init(self, radius, points, np.arange(len(payload)))
    return newinit


def _enable_pyobject_return(func):
    def getpayload(load, x):
        return (x[0], load[x[1]]) if isinstance(x, tuple) else load[x]

    @functools.wraps(func)
    def newfunc(self, *args, **kwargs):
        ret = func(self, *args, **kwargs)
        return ([getpayload(self._init_payload, r) for r in ret]
                if hasattr(ret, '__iter__')
                else getpayload(self._init_payload, ret))
    return newfunc

# when should convert point array return to tuples?
# def _enable_return_tuple(func):
#     @functools.wraps(func)
#     def newfunc(self, *args, **kwargs):
#         ret = func(self, *args, **kwargs)
#         return [(x[0],x[1],x[2]) for x in ret]
#     return newfunc


def _decorate(cls, funcname, decorator):
    orig = getattr(cls, funcname)
    setattr(cls, funcname, decorator(orig))


def _get_dtype_by_name(name):
    if name == 'object':
        return np.dtype('O')
    else:
        return vars(rif)[name].dtype


def _convert_to_arrayV3(iterable):
    seq = list(iterable)
    array = np.zeros((len(seq), 3), dtype='f4')
    for i, v in enumerate(seq):
        array[i, 0] = v[0]
        array[i, 1] = v[1]
        array[i, 2] = v[2]
    return array.view(V3)


_typemap_3d = dict()
for clsname in vars(_rif._index._stripe_index):
    cls = vars(_rif._index._stripe_index)[clsname]
    if clsname.startswith('stripe_index_3d_'):
        if clsname.endswith('_object'):
            cls.__init__ = _enable_init_pyobject_store(cls.__init__)
            _decorate(cls, '_raw_payload', _enable_pyobject_return)
            _decorate(cls, 'neighbors', _enable_pyobject_return)
        dtypes = tuple([_get_dtype_by_name(x) for x in clsname.split('_')[3:]])
        if len(dtypes) is 1:
            dtypes = (dtypes[0], None)
        _typemap_3d[dtypes] = cls
