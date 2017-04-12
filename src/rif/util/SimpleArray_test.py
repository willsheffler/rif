import numpy as np
import pytest

from rif.util import *


def test_SimpleArray_translation():
    out = echo_SimpleArray_1_int(np.ones(1))
    assert out.shape == (1,)
    assert out[0] == 1
    with pytest.raises(TypeError):
        echo_SimpleArray_1_int(np.ones(3))
    out = echo_SimpleArray_1_int(np.array([3, 2, 1])[:1])
    assert out.shape == (1,)
    assert out[0] == 3
    out = echo_SimpleArray_3_int(np.array([2, 1, 3]))
    assert out.shape == (3,)
    assert np.all(out == [2, 1, 3])
    out = echo_SimpleArray_3_int(np.array([1, 2, 3, 4, 5, 6])[1::2])
    assert out.shape == (3,)
    assert np.all(out == [2, 4, 6])


def test_SimpleArray_ref_translation():
    x = SATest()
    x.member()[0] = 3
    print(x.member())
    assert x.member()[0] == 3
