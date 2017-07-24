from rif import Atom, V3
from rif.index import stripe_index_3d
import numpy as np
import pytest


class TestAry(object):
    def __getitem__(self, i):
        return (0, 1, 2)[i]


def test_stripe_index_3d():
    points = np.arange(9, dtype='f4').view(V3)
    index = stripe_index_3d(1.0, points)
    assert np.all(index.neighbors(V3(0, 1, 2)).view('3f') == np.arange(3))
    assert np.all(index.neighbors((0, 1, 2)).view('3f') == np.arange(3))
    assert np.all(index.neighbors([0, 1, 2]).view('3f') == np.arange(3))
    assert np.all(index.neighbors(TestAry()).view('3f') == np.arange(3))
    with pytest.raises(TypeError):
        index.neighbors(object())  # no getitem
    with pytest.raises(IndexError):
        index.neighbors([1, 2])  # not enough


def test_stripe_index_3d_Atom():
    atoms = np.zeros(3, dtype=Atom)
    print(atoms['pos'])
    print(atoms)
    aindex = stripe_index_3d(3.0, atoms)


def test_stripe_index_3d_pyobject():
    points = np.arange(9, dtype='f4').view(V3)
    objs = np.array(['zero', 'one', 'two'], dtype=np.dtype(object))
    index = stripe_index_3d(1.0, points, objs)
    assert hasattr(index, '_init_payload')
    # index.__dont_gc_me = objs
    assert index.neighbor_count(V3(0, 1, 2)) == 1
    assert index.neighbor_count(V3(1, 2, 5)) == 0
    assert index.neighbors(V3(0, 1, 2))[0] is 'zero'
    assert index.neighbors(V3(3, 4, 5))[0] is 'one'
    assert index.neighbors(V3(6, 7, 8))[0] is 'two'
    assert index.neighbors((6, 7, 8))[0] is 'two'
    # no bounds check, only takes first 3 for compatibility...
    # xyzVector doesn't know how long it it
    assert index.neighbors([6, 7, 8, 'nonsense'])[0] is 'two'


def test_stripe_index_3d_payload():
    points = np.arange(9, dtype='f4').view(V3)
    atoms = np.zeros(3, dtype=Atom)
    atoms['atype'] = [1, 2, 3]
    index = stripe_index_3d(1.0, points, atoms)
    assert not hasattr(index, '_init_payload')
    # index.__dont_gc_me = objs
    assert index.neighbor_count(V3(0, 1, 2)) == 1
    assert index.neighbor_count(V3(1, 2, 5)) == 0
    # should return a numpy array????
    assert np.array_equal(index.neighbors(V3(0, 1, 2))['atype'], [1])
    assert np.array_equal(index.neighbors(V3(3, 4, 5))['atype'], [2])
    assert np.array_equal(index.neighbors(V3(6, 7, 8))['atype'], [3])
    assert np.array_equal(index.neighbors((6, 7, 8))['atype'], [3])


def test_stripe_index_3d_pyobject2():
    n = 10
    atoms = np.zeros(n, dtype=Atom)
    atoms['pos']['raw'][..., 0] = range(n)
    atoms['atype'] = range(n)
    objs = np.array([str(i) for i in range(n)], dtype=np.dtype(object))
    index = stripe_index_3d(2.5, atoms, objs)
    for i in range(len(index)):
        assert index._raw_payload(i) is objs[i]
    assert index.neighbors(atoms[0]) == list('012')
    assert index.neighbors(atoms[7]) == list('56789')


def test_stripe_index_3d_pyobject_seq_input():
    points = [(0, 1, 2), (4, 5, 6)]
    values = ['012', '456']
    index = stripe_index_3d(3.5, points, values)
    assert index.neighbors((0, 1, 2)) == ['012']
    assert index.neighbors((2, 3, 4)) == ['012', '456']
