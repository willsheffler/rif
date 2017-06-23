"""overrides numpy operators to work with rif dpypes

override numpy operators using numpy.set_numeric_ops. keep track of
overridden operators and call them if a custom operator is not found

"""
from __future__ import print_function

import sys
import numpy as np
import rif
import rif.actor
import rif.eigen_types
import pandas
from rif import V3, M3, X3, Atom
from functools import wraps

_FUNCS_BROKEN_WITH_RIOPS = [
    (np, 'cumsum'),
    (pandas.core.internals.BlockManager, 'delete')
    # (np.ndarray, 'cumsum'),  # can't set attrib of build-in ndarray!
]
_RIFOP_MODULES = [
    rif.eigen_types,
    rif.actor,
]
_NPY_RIF_OP1MAP = dict()
_NPY_RIF_OP2MAP = dict()
_ORIG_NUMPY_OPS = None

_opmap1 = dict(
    abs='absolute'
)

_opmap2 = dict(
    add='add',
    mul='multiply',
    sub='subtract',
    div='divide',
)


def _init_dispatch():
    global _NPY_RIF_OP1MAP
    global _NPY_RIF_OP2MAP
    dtmap = dict(
        fl=(type(1), type(1.0)),
        V3=V3.dtype,
        M3=M3.dtype,
        X3=X3.dtype,
        AT=Atom.dtype,
    )
    dtmap = {k: v if hasattr(v, '__iter__') else (v,)  # x -> (x,) iff not iter
             for k, v in dtmap.items()}
    for module in _RIFOP_MODULES:
        for fn in dir(module):
            if fn.startswith('rifop_'):
                splt = fn.split('_')
                if len(splt) is 3:
                    _, op, t1 = splt
                    for dt1 in dtmap[t1]:
                        k = dt1, _opmap1[op]
                        _NPY_RIF_OP1MAP[k] = getattr(module, fn)
                elif len(splt) is 4:
                    _, op, t1, t2 = splt
                    for dt1 in dtmap[t1]:
                        for dt2 in dtmap[t2]:
                            k = dt1, dt2, _opmap2[op]
                            _NPY_RIF_OP2MAP[k] = getattr(module, fn)
_init_dispatch()


def _get_type_str(t):
    try:
        return str(t.dtype)
    except AttributeError:
        return ''  # scalar


def _override1(name):
    def ufunc(x, *args, **kwargs):
        try:
            t = x.dtype if hasattr(x, 'dtype') else type(x)
            return _NPY_RIF_OP1MAP[t, name](x, *args, **kwargs)
        except (AttributeError, KeyError):
            # return getattr(np, name)(x, *args, **kwargs)
            return _ORIG_NUMPY_OPS[name](x, *args, **kwargs)
    return ufunc


def _override2(name):
    def ufunc(x, y, *args, **kwargs):
        try:
            t1 = x.dtype if hasattr(x, 'dtype') else type(x)
            t2 = y.dtype if hasattr(y, 'dtype') else type(y)
            return _NPY_RIF_OP2MAP[t1, t2, name](x, y, *args, **kwargs)
        except KeyError:
            print(name, 'x', x, 'y', y, 'args', args, 'kwargs', kwargs)
            # return getattr(np, name)(x, y, *args, **kwargs)
            return _ORIG_NUMPY_OPS[name](x, y, *args, **kwargs)
    return ufunc


def rif_operators_are_enabled():
    return _ORIG_NUMPY_OPS is not None


def global_rif_operators_enable(quiet=False):
    global _ORIG_NUMPY_OPS
    if rif_operators_are_enabled():
        print('warning: global_rif_ops is already enabled')
    else:
        d1 = {ufunc: _override1(ufunc) for ufunc in _opmap1.values()}
        d2 = {ufunc: _override2(ufunc) for ufunc in _opmap2.values()}
        d1.update(d2)
        _ORIG_NUMPY_OPS = np.set_numeric_ops(**d1)
        wrap_broken_functions()


def global_rif_operators_disable(quiet=False):
    global _ORIG_NUMPY_OPS
    assert _ORIG_NUMPY_OPS
    np.set_numeric_ops(**_ORIG_NUMPY_OPS)
    _ORIG_NUMPY_OPS = None
    unwrap_broken_functions()


class RifOperators(object):
    """contect manager for locally enabling rif ops"""

    def __enter__(self):
        self.previously_not_using_rifops = _ORIG_NUMPY_OPS is None
        if self.previously_not_using_rifops:
            global_rif_operators_enable()

    def __exit__(self, *args):
        if args:
            print("========== RifOperators: exit ============")
            for a in args:
                print(a)
            print('----------- end rifops exit --------------')
        if self.previously_not_using_rifops:
            global_rif_operators_disable()


class RifOperatorsDisabled(object):
    """context manager to locally disable RifOperators"""

    def __enter__(self):
        self.previously_using_rifops = rif_operators_are_enabled()
        if self.previously_using_rifops:
            global_rif_operators_disable()

    def __exit__(self, *args):
        if self.previously_using_rifops:
            global_rif_operators_enable()


def with_rifops_enabled(f):
    @wraps(f)
    def wrap(*args, **kwargs):
        with RifOperators():
            return f(*args, **kwargs)
    return wrap


def with_rifops_disabled(f):
    @wraps(f)
    def wrap(*args, **kwargs):
        with RifOperatorsDisabled():
            return f(*args, **kwargs)
    return wrap

_ORIG_WRAPPED_BROKEN_FUNCTIONS = None


def wrap_broken_functions():
    global _ORIG_WRAPPED_BROKEN_FUNCTIONS
    assert not _ORIG_WRAPPED_BROKEN_FUNCTIONS
    _ORIG_WRAPPED_BROKEN_FUNCTIONS = dict()
    for mod, fn in _FUNCS_BROKEN_WITH_RIOPS:
        _ORIG_WRAPPED_BROKEN_FUNCTIONS[mod.__name__, fn] = getattr(mod, fn)
        setattr(mod, fn, with_rifops_disabled(getattr(mod, fn)))


def unwrap_broken_functions():
    global _ORIG_WRAPPED_BROKEN_FUNCTIONS
    assert _ORIG_WRAPPED_BROKEN_FUNCTIONS is not None
    for mod, fn in _FUNCS_BROKEN_WITH_RIOPS:
        setattr(mod, fn, _ORIG_WRAPPED_BROKEN_FUNCTIONS[mod.__name__, fn])
    _ORIG_WRAPPED_BROKEN_FUNCTIONS = None


global_rif_operators_enable()
assert rif_operators_are_enabled()
