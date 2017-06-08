from rif import actor
from rif.actor import Atom
from rif.index import AtomIndexOneSided
import numpy as np
print(actor)


def test_AtomIndexOneSided():
    atoms = np.zeros(3, dtype=Atom)
    print(atoms['pos'])
    print(atoms)
    aindex = AtomIndexOneSided(atoms, 3.0)
    print(aindex)

    # assert 0

# NOTES:
# pyflakes doesn't validate imports
# setting python_interpreter works? at least for conda root?
# removed rif.pth from interpreters, still finding rif.actor?
# get a better handle on what's up, then post anaconda issue
