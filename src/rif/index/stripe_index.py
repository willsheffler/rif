import rif
from rif import V3, Atom, _rif
import numpy as np
import functools
import rif.dtypes
from time import clock


def stripe_index_3d(maxdis, points, payloads=None):
    """create a 3D stripe index of points[payloads] based on dtype of inputs

    stripe index is a fast index of 3D coords (possible enhanced with metadata) which can quickly get all the neighbors within a certain maxdis, which must be specified at construction

    Args:
        maxdis (float): maximum distance for points to be neighbors
        points (array): can map to array of C++ structs via dtype
        payloads (array): optional, additional data attached to each point

    Returns:
        stripe_index_3d_<points_t>_<payload_t=NULL>: Wrapper, for example, see the \*_V3 or \*_V3_V3 versions
            - :py:class:`Points only (\*_V3)
              <_rif._index._stripe_index.stripe_index_3d_V3>`
              neighbors method returns points
            - :py:class:`With payloads (\*_V3_V3)
              <_rif._index._stripe_index.stripe_index_3d_V3_V3>`
              neighbors method returns payloads (only)

    Return value notes:
        If payloads is a sequence of python objects, internally integers 1-N are
        used as the payloads, and the objects in the sequence are returned (not recommended, you should be using array style computing).

    Examples:
        - **Example: simple points with python objects**
            >>> array = np.random.uniform(size=(1000,3))
            >>> points = [(x[0], x[1], x[2]) for x in array]
            >>> type(points[13])  # converted to tuples
            <class 'tuple'>
            >>> index = stripe_index_3d(maxdis=0.1, points=points)
            >>> index
            <_rif._index._stripe_index.stripe_index_3d_V3 object at ...>
            >>> test = (0.1, 0.2, 0.3)
            >>> neighbors_fast = index.neighbor_count(test)
            >>> neighbors_slow = index.neighbor_count_brute(test)
            >>> # speedup for this example would be nothing
            >>> # because python operations dominate runtime

        - **Example: simple points with numpy**
            >>> points = np.random.uniform(size=30000).astype('f').view(V3)
            >>> points.shape
            (10000,)
            >>> index = stripe_index_3d(maxdis=0.05, points=points)
            >>> index
            <_rif._index._stripe_index.stripe_index_3d_V3 object at ...>
            >>> testpts = np.random.uniform(low=-1, high=2, size=3000).astype('f').view(V3)
            >>> tfast = clock()
            >>> neighbors_fast = index.neighbor_count(testpts)
            >>> tfast = clock() - tfast
            >>> tslow = clock()
            >>> neighbors_slow = index.neighbor_count_brute(testpts)
            >>> tslow = clock() - tslow
            >>> print('avg number of neighbors:', np.mean(neighbors_fast))  # doctest: +SKIP
            avg number of neighbors: 0.881
            >>> np.all(neighbors_fast == neighbors_slow)
            True
            >>> print('speedup is %fx vs brute force loop over points' % (tslow/tfast)) # doctest: +SKIP
            speedup is 81.904762x vs brute force loop over points
            >>> # this example is faster because it uses a numpy vectorized query
            >>> # this amortizes away the cost of python operations

        - **Example: mapping points to arbitrary python payloads**
            >>> points = [(0,0,0), (1,1,1), (2,2,2), (3,3,3)]
            >>> values = ['string payloads: ' + str(x) for x in points]
            >>> index = stripe_index_3d(maxdis=1.0, points=points, payloads=values)
            >>> index.neighbors((0,0,0))  # returns python objects
            ['string payloads: (0, 0, 0)']
            >>> index.neighbors((0.5,0.5,0.5))
            ['string payloads: (0, 0, 0)', 'string payloads: (1, 1, 1)']
            >>> # this example will also be slow... no reason to use this over brute force

        - **Example: mapping point-like structs to structs via numpy dtypes**
            >>> N = 1000
            >>> points = (np.arange(3*N).astype('f')/N).view(V3)
            >>> atoms = np.empty(N, dtype=Atom)
            >>> atoms['pos'] = points
            >>> atoms['atype'] = np.arange(N) % 23
            >>> atoms['rtype'] = np.arange(N) % 20
            >>> atoms['anum'] = np.arange(N) % 7
            >>> atoms[-3:]
            array([(([ 2.99099994,  2.9920001 ,  2.99300003],),  8, 3, 17),
                   (([ 2.99399996,  2.99499989,  2.99600005],),  9, 4, 18),
                   (([ 2.99699998,  2.99799991,  2.99900007],), 10, 5, 19)],
                  dtype=[('pos', [('raw', '<f4', (3,))]), ('atype', 'u1'), ('anum', 'u1'), ('rtype', '<u2')])
            >>> index = stripe_index_3d(maxdis=0.05, points=points, payloads=atoms)
            >>> index.neighbors((1, 1, 1))
            array([(([ 0.972     ,  0.97299999,  0.97399998],),  2, 2,  4),
                   (([ 0.97500002,  0.97600001,  0.977     ],),  3, 3,  5),
                   (([ 0.97799999,  0.97899997,  0.98000002],),  4, 4,  6),
                   (([ 0.98100001,  0.98199999,  0.98299998],),  5, 5,  7),
                   (([ 0.98400003,  0.98500001,  0.986     ],),  6, 6,  8),
                   (([ 0.98699999,  0.98799998,  0.98900002],),  7, 0,  9),
                   (([ 0.99000001,  0.991     ,  0.99199998],),  8, 1, 10),
                   (([ 0.99299997,  0.99400002,  0.995     ],),  9, 2, 11),
                   (([ 0.99599999,  0.99699998,  0.99800003],), 10, 3, 12),
                   (([ 0.99900001,  1.        ,  1.00100005],), 11, 4, 13),
                   (([ 1.00199997,  1.00300002,  1.00399995],), 12, 5, 14),
                   (([ 1.005     ,  1.00600004,  1.00699997],), 13, 6, 15),
                   (([ 1.00800002,  1.00899994,  1.00999999],), 14, 0, 16),
                   (([ 1.01100004,  1.01199996,  1.01300001],), 15, 1, 17),
                   (([ 1.01400006,  1.01499999,  1.01600003],), 16, 2, 18),
                   (([ 1.01699996,  1.01800001,  1.01900005],), 17, 3, 19),
                   (([ 1.01999998,  1.02100003,  1.02199996],), 18, 4,  0),
                   (([ 1.023     ,  1.02400005,  1.02499998],), 19, 5,  1),
                   (([ 1.02600002,  1.02699995,  1.028     ],), 20, 6,  2)],
                  dtype=[('pos', [('raw', '<f4', (3,))]), ('atype', 'u1'), ('anum', 'u1'), ('rtype', '<u2')])
    """
    dtype2 = None
    if payloads is not None:
        assert len(payloads) == len(points)
        dtype2 = np.dtype('O')
        if hasattr(payloads, 'dtype'):
            dtype2 = payloads.dtype

    dtype1 = np.dtype('O')
    if hasattr(points, 'dtype'):
        dtype1 = points.dtype
    else:
        # if payloads is sequence of objects, allow points
        # to be an iterable of 'index-ables' of len at least 3
        # for example, a list of rosetta xyzVectors
        points = _convert_to_arrayV3(points)
        dtype1 = V3.dtype
    try:
        return _typemap_3d[dtype1, dtype2](maxdis, points, payloads)
    except KeyError:
        raise ValueError(
            "stripe_index_3d only supports numpy arrays (can have " +
            "fancy dtypes) or object payloads (which is inefficient), not: " +
            str(type(points)) + ', ' + str(type(payloads)))


def _enable_init_pyobject_store(init):
    @functools.wraps(init)
    def newinit(self, maxdis, points, payloads):
        self._init_payload = list(payloads)
        init(self, maxdis, points, np.arange(len(payloads)))
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
