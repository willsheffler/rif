import numpy as np

from rif.dtypes import rif_ops
from rif.eigen_types import v3f_t
from rif.actor import atom_t


def test_atom_dtype():
    a = np.zeros(10, dtype=atom_t).view(np.recarray)
    assert a.shape == (10,)
    assert a.pos['crd'].shape == (10, 3)
    assert a.atype.shape == (10,)
    assert a.rtype.shape == (10,)
    assert a.anum.shape == (10,)
    a.pos['crd'] = np.random.randn(10, 3)


def test_atom_math():
    v = np.ones(2, dtype=v3f_t)
    a = np.zeros(2, dtype=atom_t)
    a['atype'] = [3, 4]
    a['rtype'] = [12, 6]
    a['anum'] = [1, 2]
    print('v', v)
    print('a', a)
    with rif_ops():
        assert np.all(1 == (a + v)['pos']['crd'])
        assert np.all(3 == (v + a + v + v)['pos']['crd'])
        assert np.all(-1 == (a - v)['pos']['crd'])
        assert np.all(3 == (3 * v - a + v - v)['pos']['crd'])
