from rif.numeric.lattice import *
from rif.util import sa_dtype
import numpy as np
import pytest


def test_bcc():
    N = 3
    ncell = np.array(N * [4])
    lower = np.array(N * [-10])
    upper = np.array(N * [10])
    bcc = BCC3(nc=ncell, lb=lower, ub=upper)
    idx = np.arange(2)
    bincen = bcc.center(idx)
    assert bincen['raw'].shape == (2, N)
    # print(bincen)
    crd = np.ones((2, N))
    assert bcc.index(crd).shape == (2,)
    assert bcc.index(crd.view(sa_dtype[3, 'f8'])).shape == (2,)


def test_bcc_invertability():
    N = 6
    bcc = BCC6(N * [20], N * [-20], N * [20])
    idx = np.arange(1000)
    bincen = bcc.center(idx)
    idx2 = bcc.index(bincen['raw'])
    assert np.all(idx == idx2)


def test_bcc_bounds():
    for nc0 in range(3, 10):
        qsb = +1.0 + 2.0 / (nc0 - 2.0)
        nc = np.array([3, 3, 3, nc0, nc0])
        ub = np.array([10, 10, 10, qsb, qsb])
        bcc = BCC5(nc=nc, lb=-ub, ub=ub)
        s = set()
        for i in range(bcc.ncells):
            s.add(bcc.center(i)[3])
            s.add(bcc.center(i)[4])
        # print(nc0, qsb, sorted(s))
        s = np.array(sorted(s))
        print(nc0, s)
        assert np.min(np.abs(s - 1)) < 0.000001  # have 1.0
        assert np.min(np.abs(s + 1)) < 0.000001  # have -1.0

        assert len(s) == 2 * nc[-1]
