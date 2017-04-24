from __future__ import print_function
import pytest
from rif.data.rotamers import richardson_rots as rr
from rif.util import rcl
import numpy as np

# @pytest.mark.skip('not implemented')
# @pytest.mark.skipif('not rcl.HAVE_PYROSETTA')


@pytest.mark.skip(reason='richardson rots not implemented yet')
def test_rotamers():
    print(rr)
    assert rr.shape[0] == 7
    assert rr.shape[1] == 9
    print(rr.chi1[4])
    print("=== just listed rotamers ===")
    for i in range(rr.shape[0]):
        print('chi1/2', rr.chi1[i], rr.chi2[i])
    print("=== chi2 every 20 degrees ===")
    for i in range(rr.shape[0]):
        for chi2 in np.arange(rr.lb2[i], rr.ub2[i], 20.0):
            print('chi1/2', rr.chi1[i], chi2)
    # assert 0
