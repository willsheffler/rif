#include <iostream>
#include "stdint.h"

#include "numpy/arrayobject.h"
#include "pybind11/numpy.h"
#include "pybind11/pybind11.h"

#include "ndarray.h"
#include "ndarray_pybind11_converter.hh"

struct ExampleStruct {
  double a;
  int b;
};

struct NestedExampleStruct {
  ExampleStruct ex1;
  ExampleStruct ex2;
  bool flag;
};

void sum_example_struct(ndarray::Array<ExampleStruct, 1> input,
                        ndarray::Array<double, 1> out) {
  if (input.getSize<0>() != out.getSize<0>()) {
    throw std::invalid_argument("Invalid array shapes");
  }

  for (std::size_t i = 0; i < input.getSize<0>(); ++i) {
    out[i] = input[i].a + input[i].b;
  }
}

void sum_example(ndarray::Array<double, 1> a, ndarray::Array<double, 1> b,
                 ndarray::Array<double, 1> out) {
  if ((a.getSize<0>() != b.getSize<0>()) ||
      (a.getSize<0>() != out.getSize<0>())) {
    throw std::invalid_argument("Invalid array shapes");
  }

  for (std::size_t i = 0; i < a.getSize<0>(); ++i) {
    out[i] = a[i] + b[i];
  }
}

namespace py = pybind11;

void RIFLIB_PYBIND_test_ndarray(py::module &m) {
  // _import_array initializes numpy api, needed for ndarray/array
  // interconversion

  if (_import_array() < 0) {
    PyErr_SetString(PyExc_ImportError,
                    "numpy.core.multiarray failed to import");
    std::cerr << "numpy.core.multiarray failed to import" << std::endl;
    return;
  };

  PYBIND11_NUMPY_DTYPE(ExampleStruct, a, b);

  m.def("sum_example", &sum_example);
  m.def("sum_example_struct", &sum_example_struct);
  m.attr("real_dtype") = py::dtype::of<double>();
  m.attr("struct_dtype") = py::dtype::of<ExampleStruct>();
}
