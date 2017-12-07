import rif.dtypes

from _rif._test import _pandas_interop as pi

import numpy as np
import xarray as xr


# def test_can_treat_xr_as_dict():
    # data = xr.DataArray(np.random.randn(2, 3), coords={
        # 'x': ['a', 'b']}, dims=('x', 'y'))
    # ds = xr.Dataset({'foo': data, 'bar': ('x', [1, 2]), 'baz': np.pi})
    # sp = pi.print_dict(ds)
    # print(sp)
    # assert 17 == len(sp.splitlines())
    # kn, totlen = pi.key_names_and_len(ds)
    # assert totlen == 11


class C(np.ndarray):
    "needs to store numpy.index_exp[1:3, 1:3] and array to set this to"
    def __new__(cls, *args, **kwargs):
        print('In __new__ with class %s' % cls)
        return super(C, cls).__new__(cls, *args, **kwargs)

    def __init__(self, *args, **kwargs):
        # in practice you probably will not need or want an __init__
        # method for your subclass
        print('In __init__ with class %s' % self.__class__)

    def __array_finalize__(self, obj):
        print('In __array_finalize____:')
        print('   self type is %s' % type(self))
        print('   obj type is %s' % type(obj))
        if obj is None:
            return

    def __array_ufunc__(self, ufunc, method, *inputs, **kwargs):
        args, in_no = [], []
        for i, input_ in enumerate(inputs):
            if isinstance(input_, A):
                in_no.append(i)
                args.append(input_.view(np.ndarray))
            else:
                args.append(input_)

        outputs = kwargs.pop('out', None)
        out_no = []
        if outputs:
            out_args = []
            for j, output in enumerate(outputs):
                if isinstance(output, A):
                    out_no.append(j)
                    out_args.append(output.view(np.ndarray))
                else:
                    out_args.append(output)
            kwargs['out'] = tuple(out_args)
        else:
            outputs = (None,) * ufunc.nout

        info = {}
        if in_no:
            info['inputs'] = in_no
        if out_no:
            info['outputs'] = out_no

        results = super(A, self).__array_ufunc__(ufunc, method,
                                                 *args, **kwargs)
        if results is NotImplemented:
            return NotImplemented

        if method == 'at':
            if isinstance(inputs[0], A):
                inputs[0].info = info
            return

        if ufunc.nout == 1:
            results = (results,)

        results = tuple((np.asarray(result).view(A)
                         if output is None else output)
                        for result, output in zip(results, outputs))
        if results and isinstance(results[0], A):
            results[0].info = info

        return results[0] if len(results) == 1 else results


def test_ndarray_subclass():
    a = np.zeros((3,))
    ca = a.view(C)
    assert type(ca) is C
    v = ca[1:]
    assert type(v) is C
    # assert 0


def test_xarray_with_ndarray_subclass():
    sub = np.random.randn(2, 3).view(C)
    assert type(sub) is C
    # x = xr.DataArray(
    # sub,
    # coords={'x': ['a', 'b'], 'y': ['c', 'd']},
    # name='somename',
    # fastpath=True)
    # assert type(x._variable) is C
    # assert type(x.data) is C
    # assert 0
