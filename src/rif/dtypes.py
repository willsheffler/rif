from __future__ import print_function

import sys
import numpy as np
import rif.eigen_types
from rif_cpp.dtypes import print_numpy_info
from rif.eigen_types import *
from rif.actor import Atom


def init_dispatch():
    tM = str(m3f_t)
    tV = str(v3f_t)
    tA = str(Atom.dtype)

    # str is necessary for python2
    # TODO: figures out why my dtypes won't hash in python2
    npy_rif_op1map = dict()
    npy_rif_op1map[tM, 'absolute'] = rif.eigen_types.abs_m3f
    npy_rif_op1map[tV, 'absolute'] = rif.eigen_types.abs_v3f

    npy_rif_op2map = dict()
    npy_rif_op2map[tV, tV, 'add'] = rif.eigen_types.add_v3f
    npy_rif_op2map[tV, tV, 'subtract'] = rif.eigen_types.sub_v3f
    npy_rif_op2map[tV, tV, 'multiply'] = rif.eigen_types.mul_v3f
    npy_rif_op2map[tM, tM, 'add'] = rif.eigen_types.add_m3f
    npy_rif_op2map[tV, '', 'multiply'] = rif.eigen_types.mul_v3f_f
    npy_rif_op2map[tV, '', 'divide'] = rif.eigen_types.div_v3f_f
    npy_rif_op2map['', tV, 'multiply'] = rif.eigen_types.mul_f_v3f
    npy_rif_op2map[tM, '', 'multiply'] = rif.eigen_types.mul_m3f_f
    npy_rif_op2map[tM, '', 'divide'] = rif.eigen_types.div_m3f_f
    npy_rif_op2map['', tM, 'multiply'] = rif.eigen_types.mul_f_m3f
    npy_rif_op2map[tM, tM, 'subtract'] = rif.eigen_types.sub_m3f
    npy_rif_op2map[tM, tM, 'multiply'] = rif.eigen_types.mul_m3f
    npy_rif_op2map[tM, tV, 'multiply'] = rif.eigen_types.mul_m3f_v3f
    npy_rif_op2map[tA, tV, 'add'] = rif.actor.add_atom_v3f
    npy_rif_op2map[tV, tA, 'add'] = rif.actor.add_v3f_atom
    npy_rif_op2map[tA, tV, 'subtract'] = rif.actor.sub_atom_v3f
    npy_rif_op2map[tV, tA, 'subtract'] = rif.actor.sub_v3f_atom
    return npy_rif_op1map, npy_rif_op2map


npy_rif_op1map, npy_rif_op2map = init_dispatch()


def get_type_str(t):
    try:
        return str(t.dtype)
    except AttributeError:
        return ''  # scalar


def override1(name):
    def ufunc(x):
        try:
            t = get_type_str(x)
            r = npy_rif_op1map[t, name](x)
            return r
        except (AttributeError, KeyError):
            r = getattr(np, name)(x)
            return r
    return ufunc


def override2(name):
    def ufunc(x, y):
        try:
            t1 = get_type_str(x)
            t2 = get_type_str(y)
            r = npy_rif_op2map[t1, t2, name](x, y)
            return r
        except (AttributeError, KeyError):
            r = getattr(np, name)(x, y)
            return r
    return ufunc


_GLOBAL_RIF_OPS = None


class rif_ops(object):
    """contect manager for locally enabling rif ops"""

    def __enter__(self):
        global _GLOBAL_RIF_OPS
        if _GLOBAL_RIF_OPS is not None:
            print('warning: global_rif_ops is enabled, rif_ops contect manager doing nothing')
            self.orig = None
        else:
            print('rif_ops: enter')
            d = {ufunc: override1(ufunc) for ufunc in ('absolute'.split())}
            d2 = {ufunc: override2(ufunc)
                  for ufunc in ('add subtract multiply'.split())}
            d.update(d2)
            self.orig = np.set_numeric_ops(**d)

    def __exit__(self, *args):
        if args:
            print("rif_ops: exit", args)
        if self.orig:
            np.set_numeric_ops(**self.orig)


class rif_ops_disable(object):
    """context manager to locally disable rif_ops"""

    def __enter__(self):
        global_rif_ops_disable(quiet=True)

    def __exit__(self, *args):
        global_rif_ops_enable(quiet=True)


def global_rif_ops_enable(quiet=False):
    global _GLOBAL_RIF_OPS
    if not _GLOBAL_RIF_OPS:
        tmp = rif_ops()
        tmp.__enter__()
        _GLOBAL_RIF_OPS = tmp
    elif not quiet:
        print("warning: global_rif_ops_enable called after already enabled")


def global_rif_ops_disable(quiet=False):
    global _GLOBAL_RIF_OPS
    if _GLOBAL_RIF_OPS:
        tmp = _GLOBAL_RIF_OPS
        _GLOBAL_RIF_OPS = None
        tmp.__exit__()
    elif not quiet:
        print("warning: global_rif_ops_disable called when not enabled")


if sys.version_info.major is 3:
    global_rif_ops_enable()
