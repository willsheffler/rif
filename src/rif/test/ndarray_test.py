from rif_cpp.test import ndarray

import numpy
import numpy.testing


def test_basic():
    inp = numpy.arange(10, dtype=numpy.float)
    out = numpy.empty(10)

    ndarray.sum_example(inp, inp, out)
    numpy.testing.assert_almost_equal(inp + inp, out)


def test_struct():
    inp = numpy.empty(10, dtype=ndarray.struct_dtype)
    inp["a"] = numpy.arange(10)
    inp["b"] = 10

    out = numpy.empty(10)

    ndarray.sum_example_struct(inp, out)
    numpy.testing.assert_almost_equal(inp["a"] + inp["b"], out)
