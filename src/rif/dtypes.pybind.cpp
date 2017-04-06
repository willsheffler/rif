#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

#include <iostream>

namespace py = pybind11;

void print_numpy_info(py::array a) {
  py::buffer_info b(py::buffer(a).request());
  std::cout << "size: " << b.size << std::endl;
  std::cout << "itemsize: " << b.itemsize << std::endl;
  std::cout << "ndim: " << b.ndim << std::endl;
  std::cout << "format: '" << b.format << "'" << std::endl;
}

void RIFLIB_PYBIND_rif_dtypes(py::module& m) {
  m.def("print_numpy_info", &print_numpy_info);
}
