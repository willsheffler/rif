#include "pybind11/numpy.h"
#include "pybind11/pybind11.h"

#include "eigen_types.hpp"

namespace py = pybind11;
using namespace rif;

namespace std {
template <class F, int M, int N>
struct is_pod<Eigen::Matrix<F, M, N>>
    : public std::integral_constant<bool, true> {};
}

V3f test_add_v3f(V3f a, V3f b) { return a + b; }

void RIFLIB_PYBIND_eigen_types(py::module &m) {
  ::pybind11::detail::npy_format_descriptor<V3f>::register_dtype(
      {::pybind11::detail::field_descriptor{
           "x", 0, sizeof(float), alignof(float),
           ::pybind11::format_descriptor<float>::format(),
           ::pybind11::detail::npy_format_descriptor<float>::dtype()},
       ::pybind11::detail::field_descriptor{
           "y", 4, sizeof(float), alignof(float),
           ::pybind11::format_descriptor<float>::format(),
           ::pybind11::detail::npy_format_descriptor<float>::dtype()},
       ::pybind11::detail::field_descriptor{
           "z", 8, sizeof(float), alignof(float),
           ::pybind11::format_descriptor<float>::format(),
           ::pybind11::detail::npy_format_descriptor<float>::dtype()}});

  m.attr("v3f_t") = py::dtype::of<V3f>();
  m.def("add_v3f_test", py::vectorize(test_add_v3f));
}
