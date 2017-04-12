from rif.numeric.bcc_lattice import *
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
    assert bincen.shape == (2, N)
    # print(bincen)
    crd = np.ones((2, N))
    assert bcc.index(crd).shape == (2,)
    assert bcc.index(crd.view(sa_dtype[3, 'f8'])).shape == (2,)


def test_bcc_invertability():
    N = 6
    bcc = BCC6(N * [20], N * [-20], N * [20])
    idx = np.arange(1000)
    bincen = bcc.center(idx)
    idx2 = bcc.index(bincen)
    assert np.all(idx == idx2)
