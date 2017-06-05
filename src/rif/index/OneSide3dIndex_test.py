from rif.actor import Atom
from rif.index import AtomIndexOneSided
import numpy as np


def test_AtomIndexOneSided():
    atoms = np.zeros(3, dtype=Atom)
    print(atoms['pos'])
    print(atoms)
    aindex = AtomIndexOneSided(atoms, 3.0)
    print(aindex)

    # assert 0
