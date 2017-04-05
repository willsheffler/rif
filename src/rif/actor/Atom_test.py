import numpy as np

from rif.actor import atom_t


def test_atom_dtype():
    a = np.zeros(10, dtype=atom_t).view(np.recarray)
    a.x = np.random.randn(10)
    assert a.shape == (10,)
    assert a.x.shape == (10,)
    assert a.y.shape == (10,)
    assert a.z.shape == (10,)
    assert a.atype.shape == (10,)
    assert a.rtype.shape == (10,)
    assert a.anum.shape == (10,)
