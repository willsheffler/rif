import pytest
import os
from os.path import join, dirname, abspath, exists


@pytest.fixture(scope='session')
def rootdir():
    d = abspath(dirname(__file__))
    assert exists(d)
    return d


@pytest.fixture(scope='session')
def datadir(rootdir):
    d = join(rootdir, 'data')
    assert exists(d)
    return d


@pytest.fixture(scope='session')
def pdbfname(datadir):
    f = join(datadir, 'pdb', '3uc7A.pdb')
    assert exists(f)
    return f
