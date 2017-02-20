#include <pybind11/pybind11.h>
#include <iostream>

void np_array_info() { std::cout << "np_array_info" << std::endl; }

namespace py = pybind11;

void RIFLIB_PYBIND_test_numpy(py::module &m) {
  m.def("np_array_info", &np_array_info);
}
