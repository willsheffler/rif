import numpy as np
import rif
from rif.eigen_types import V3


# does this fine need to exist???

def test_np_info():
    x = np.ones(5, dtype=[('raw', '(2,3)f4')])
    rif.eigen_types.print_numpy_info(x)


def test_V3():
    assert V3(1, 2, 3) == V3([1, 2, 3])
    assert V3(1, 2, 3) != V3([1, 2, 4])
