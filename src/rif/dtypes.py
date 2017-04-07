from __future__ import print_function

import sys
import numpy as np
import rif.eigen_types
from rif_cpp.dtypes import print_numpy_info
from rif.eigen_types import v3f_t, m3f_t, v3i_t
from rif.actor import atom_t

# str is necessary for python2
# TODO: figures out why my dtypes won't hash in python2
npy_rif_op1map = dict()
npy_rif_op1map[str(m3f_t), 'absolute'] = rif.eigen_types.abs_m3f
npy_rif_op1map[str(v3f_t), 'absolute'] = rif.eigen_types.abs_v3f

npy_rif_op2map = dict()
npy_rif_op2map[str(v3f_t), str(v3f_t), 'add'] = rif.eigen_types.add_v3f
npy_rif_op2map[str(v3f_t), str(v3f_t), 'subtract'] = rif.eigen_types.sub_v3f
npy_rif_op2map[str(v3f_t), str(v3f_t), 'multiply'] = rif.eigen_types.mul_v3f
npy_rif_op2map[str(m3f_t), str(m3f_t), 'add'] = rif.eigen_types.add_m3f
npy_rif_op2map[str(v3f_t), '', 'multiply'] = rif.eigen_types.mul_v3f_f
npy_rif_op2map[str(v3f_t), '', 'divide'] = rif.eigen_types.div_v3f_f
npy_rif_op2map['', str(v3f_t), 'multiply'] = rif.eigen_types.mul_f_v3f
npy_rif_op2map[str(m3f_t), '', 'multiply'] = rif.eigen_types.mul_m3f_f
npy_rif_op2map[str(m3f_t), '', 'divide'] = rif.eigen_types.div_m3f_f
npy_rif_op2map['', str(m3f_t), 'multiply'] = rif.eigen_types.mul_f_m3f
npy_rif_op2map[str(m3f_t), str(m3f_t), 'subtract'] = rif.eigen_types.sub_m3f
npy_rif_op2map[str(m3f_t), str(m3f_t), 'multiply'] = rif.eigen_types.mul_m3f
npy_rif_op2map[str(m3f_t), str(v3f_t), 'multiply'] = rif.eigen_types.mul_m3f_v3f
npy_rif_op2map[str(atom_t), str(v3f_t), 'add'] = rif.actor.add_atom_v3f
npy_rif_op2map[str(v3f_t), str(atom_t), 'add'] = rif.actor.add_v3f_atom
npy_rif_op2map[str(atom_t), str(v3f_t), 'subtract'] = rif.actor.sub_atom_v3f
npy_rif_op2map[str(v3f_t), str(atom_t), 'subtract'] = rif.actor.sub_v3f_atom


def get_type(t):
    try:
        return str(t.dtype)
    except AttributeError:
        return ''  # scalar


def override1(name):
    def ufunc(x):
        try:
            t = get_type(x)
            r = npy_rif_op1map[t, name](x)
            return r
        except (AttributeError, KeyError):
            r = getattr(np, name)(x)
            return r
    return ufunc


def override2(name):
    def ufunc(x, y):
        try:
            t1 = get_type(x)
            t2 = get_type(y)
            r = npy_rif_op2map[t1, t2, name](x, y)
            return r
        except (AttributeError, KeyError):
            r = getattr(np, name)(x, y)
            return r
    return ufunc


class rif_ops(object):
    def __enter__(self):
        print('rif_ops: enter')
        d = {ufunc: override1(ufunc) for ufunc in ('absolute'.split())}
        d2 = {ufunc: override2(ufunc) for ufunc in ('add subtract multiply'.split())}
        d.update(d2)
        self.orig = np.set_numeric_ops(**d)

    def __exit__(self, *args):
        print("rif_ops: exit", args)
        np.set_numeric_ops(**self.orig)
