from rif.numeric.bcc_lattice import *
import numpy as np


def test_bcc():
    bcc = BCC10f4u8(ncell=np.array(10 * [4]),
                    lb=np.array(10 * [-10]),
                    ub=np.array(10 * [10]))
    # assert 0
