#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <iostream>
#include <sstream>

#include "dict_testing.hpp"

namespace py = pybind11;

void RIFLIB_PYBIND_test_xarray_interop(py::module &m) {
  m.def("print_dict", &rif::test::print_dict);
  m.def("key_names_and_len", &rif::test::key_names_and_len);
}
