#include "pybind11/numpy.h"
#include "pybind11/pybind11.h"

namespace py = pybind11;

struct STRUCT {
  float FIELD1;
  float FIELD2;
};

void RIFLIB_PYBIND_test_dtype_test(py::module &m) {
  // PYBIND11_NUMPY_DTYPE(STRUCT, FIELD1);
  ::pybind11::detail::npy_format_descriptor<STRUCT>::register_dtype(
      {::pybind11::detail::field_descriptor{
           "FIELD1", 0, sizeof(float), alignof(float),
           ::pybind11::format_descriptor<float>::format(),
           ::pybind11::detail::npy_format_descriptor<float>::dtype()},
       ::pybind11::detail::field_descriptor{
           "FIELD2", 4, sizeof(float), alignof(float),
           ::pybind11::format_descriptor<float>::format(),
           ::pybind11::detail::npy_format_descriptor<float>::dtype()}});

  m.attr("struct_t") = py::dtype::of<STRUCT>();
}
