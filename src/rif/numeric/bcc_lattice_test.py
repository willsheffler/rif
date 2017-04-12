from rif.numeric.bcc_lattice import *
import numpy as np
import pytest


def test_bcc():
    bcc = BCC10f4u8(ncell=np.array(10 * [4]),
                    lb=np.array(10 * [-10]),
                    ub=np.array(10 * [10]))
    print(bcc[np.arange(2)].shape)
    with pytest.raises(ValueError):
        bcc = BCC10f4u4(ncell=np.array(10 * [1000]),
                        lb=np.array(10 * [-10]),
                        ub=np.array(10 * [10]))
