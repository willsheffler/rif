from __future__ import print_function

import sys
import numpy as np
import rif
from rif import V3, M3, X3, Atom
import rif.actor


def _init_dispatch():
    _npy_rif_op1map = dict()
    _npy_rif_op2map = dict()
    if sys.version_info.major > 2:
        modules_to_search = [rif.eigen_types, rif.actor]

        _npy_rif_op1map[M3.dtype, 'absolute'] = rif.eigen_types.op_abs_M3
        _npy_rif_op1map[V3.dtype, 'absolute'] = rif.eigen_types.op_abs_V3

        # binary ops set here

        opmap = dict(
            add='add',
            mul='multiply',
            sub='subtract',
            div='divide')
        dtmap = dict(
            fl=(type(1), type(1.0)),
            V3=V3.dtype,
            M3=M3.dtype,
            X3=X3.dtype,
            AT=Atom.dtype,
        )
        dtmap = {k: v if hasattr(v, '__iter__') else (v,)
                 for k, v in dtmap.items()}

        for module in modules_to_search:
            for fn in dir(module):
                splt = fn.split('_')
                if len(splt) is not 4:
                    continue
                head, op, t1, t2 = splt
                if head != "op":
                    continue
                for dt1 in dtmap[t1]:
                    for dt2 in dtmap[t2]:
                        k = dt1, dt2, opmap[op]
                        # print(k, fn)
                        _npy_rif_op2map[k] = getattr(module, fn)
    return _npy_rif_op1map, _npy_rif_op2map


_npy_rif_op1map, _npy_rif_op2map = _init_dispatch()


def _get_type_str(t):
    try:
        return str(t.dtype)
    except AttributeError:
        return ''  # scalar


def _override1(name):
    def ufunc(x):
        try:
            t = x.dtype if hasattr(x, 'dtype') else type(x)
            r = _npy_rif_op1map[t, name](x)
            return r
        except (AttributeError, KeyError):
            r = getattr(np, name)(x)
            return r
    return ufunc


def _override2(name):
    def ufunc(x, y):
        try:
            t1 = x.dtype if hasattr(x, 'dtype') else type(x)
            t2 = y.dtype if hasattr(y, 'dtype') else type(y)
            r = _npy_rif_op2map[t1, t2, name](x, y)
            return r
        except KeyError:
            r = getattr(np, name)(x, y)
            return r
    return ufunc


_GLOBAL_RIF_OPS = None


class RifOperators(object):
    """contect manager for locally enabling rif ops"""

    def __enter__(self):
        if sys.version_info.major == 2:
            print("RifOperators not supported in python 2")
            return
        global _GLOBAL_RIF_OPS
        if _GLOBAL_RIF_OPS is not None:
            print('warning: global_rif_ops is enabled')
            print('RifOperators context manager doing nothing')
            self.orig = None
        else:
            print('RifOperators: enter')
            d = {ufunc: _override1(ufunc) for ufunc in ('absolute'.split())}
            d2 = {ufunc: _override2(ufunc)
                  for ufunc in ('add subtract multiply divide'.split())}
            d.update(d2)
            self.orig = np.set_numeric_ops(**d)

    def __exit__(self, *args):
        if sys.version_info.major == 2:
            print("RifOperators not supported in python 2")
            return

        if args:
            print("RifOperators: exit", args)
        if self.orig:
            np.set_numeric_ops(**self.orig)


class RifOperatorsDisabled(object):
    """context manager to locally disable RifOperators"""

    def __enter__(self):
        global_rif_operators_disable(quiet=True)

    def __exit__(self, *args):
        global_rif_operators_enable(quiet=True)


def global_rif_operators_enable(quiet=False):
    if sys.version_info.major == 2:
        print("RifOperators not supported in python 2")
        return
    global _GLOBAL_RIF_OPS
    if not _GLOBAL_RIF_OPS:
        tmp = RifOperators()
        tmp.__enter__()
        _GLOBAL_RIF_OPS = tmp
    elif not quiet:
        print("warning: global_rif_operators_enable called after already enabled")


def global_rif_operators_disable(quiet=False):
    if sys.version_info.major == 2:
        print("RifOperators not supported in python 2")
        return
    global _GLOBAL_RIF_OPS
    if _GLOBAL_RIF_OPS:
        tmp = _GLOBAL_RIF_OPS
        _GLOBAL_RIF_OPS = None
        tmp.__exit__()
    elif not quiet:
        print("warning: global_rif_operators_disable called when not enabled")


if sys.version_info.major is 3:
    global_rif_operators_enable()
