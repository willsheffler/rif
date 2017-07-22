from rif import Atom, V3
from rif.index import stripe_index_3d
import numpy as np
import pytest


def test_stripe_index_3d_Atom():
    atoms = np.zeros(3, dtype=Atom)
    print(atoms['pos'])
    print(atoms)
    aindex = stripe_index_3d(3.0, atoms)

    print(aindex)


@pytest.mark.xfail
def test_ObjectIndexOneSided():
    points = np.arange(9, dtype='f4').view(V3)
    objs = np.array(['zero', 'one', 'two'], dtype=np.dtype(object))
    index = stripe_index_3d(1.0, points, objs)
    assert hasattr(index, '_init_payload')
    # index.__dont_gc_me = objs
    assert index.nbcount(V3(0, 1, 2)) == 1
    assert index.nbcount(V3(1, 2, 5)) == 0
    # for i in range(len(index)):
    # print('test_get:', index._raw(i))
    # del index
    # print(objs)
    assert 0


@pytest.mark.xfail
def test_Atom2ObjectIndexOneSided():
    atoms = np.zeros(3, dtype=Atom)
    atoms['atype'] = [1, 2, 3]
    print(atoms)
    objs = np.array(['zero', 'one', 'two'], dtype=np.dtype(object))
    # for o in objs:
    # print(id(o))
    index = Atom2ObjectIndexOneSided(3.0, atoms, objs)

    for i in range(len(index)):
        print('test_get:', index._raw_payload(i), id(index._raw_payload(i)))
        assert index._raw_payload(i) is objs[i]


# NOTES:
# pyflakes doesn't validate imports
# setting python_interpreter works? at least for conda root?
# removed rif.pth from interpreters, still finding rif.actor?
# get a better handle on what's up, then post anaconda issue
