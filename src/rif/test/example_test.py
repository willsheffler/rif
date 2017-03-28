import rif


def test_pybind_example():
    print(rif.__file__)
    assert rif.test.example.add(1, 2) == 3
    assert rif.test.example.sub(1, 2) == -1
    assert rif.test.example.mul(1, 2) == 2
