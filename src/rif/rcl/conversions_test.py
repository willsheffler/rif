from rif.rcl.conversions import *
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


@pytest.mark.skipif('not rcl.HAVE_PYROSETTA')
def test_rays(pose):
    a = rays(pose, 'n->h')
    assert a.shape == (6,)  # 1st res fails Nterm 1H not H
    a = a.view('(4,2)f')
    d = a[..., :3, 1]
    d = (d * d).sum(axis=1)
    assert np.all(np.abs(d - 1.0) < 0.001)

    b = rays(pose, 'c->o')
    assert b.shape == (7,)

    double = rays(pose, 'c->o n->h', shifts=[0, 1])
    assert double.shape == (6, 2)


@pytest.mark.skipif('not rcl.HAVE_PYROSETTA')
def test_rosetta_stub_to_matrix():
    stub = rcl.Stub()
    stub.M.xy = 1
    stub.v.x = 10
    stub.v.y = 11
    stub.v.z = 12
    rifstub = to_rif_stub(stub)
    stub144 = rifstub.view('4,4f')
    assert np.all(stub144[0, 0, :] == [1, 1, 0, 10])
    assert np.all(stub144[0, 1, :] == [0, 1, 0, 11])
    assert np.all(stub144[0, 2, :] == [0, 0, 1, 12])
    assert np.all(stub144[0, 3, :] == [0, 0, 0, 1])
    rosstub = to_rosetta_stub(rifstub)
    assert rosstub == stub
