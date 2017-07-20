from rif.rcl.conversions import atoms, stubs, rays
from rif.chem.biochem import check_atom_types
import numpy as np
from rif import rcl
import pytest


@pytest.mark.skipif('not rcl.HAVE_PYROSETTA')
def test_atoms_fixed_width(pose):
    a = atoms(pose, 'n ca c o (cb or ca)')
    assert a.shape == (7, 5)
    check_atom_types(a)


@pytest.mark.skipif('not rcl.HAVE_PYROSETTA')
def test_atoms_all(pose):
    a = atoms(pose, 'all')
    assert a.shape == (104,)
    check_atom_types(a)


@pytest.mark.skipif('not rcl.HAVE_PYROSETTA')
def test_atoms_heavy(pose):
    a = atoms(pose, 'heavy')
    assert a.shape == (58,)
    check_atom_types(a)
