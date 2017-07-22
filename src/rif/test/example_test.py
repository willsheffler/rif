import _rif


def test_pybind_example():
    print(_rif.__file__)
    assert _rif._test._example.add(1, 2) == 3
    assert _rif._test._example.sub(1, 2) == -1
    assert _rif._test._example.mul(1, 2) == 2
