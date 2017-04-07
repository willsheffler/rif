#include "util/pybind_numpy.hpp"

#include "eigen_types.hpp"

#include <iostream>

namespace py = pybind11;
using namespace rif;

V3f add_v3f(V3f a, V3f b) { return a + b; }
V3f sub_v3f(V3f a, V3f b) { return a - b; }
float mul_v3f(V3f a, V3f b) { return a.dot(b); }
float abs_v3f(V3f a) { return a.norm(); }
M3f add_m3f(M3f a, M3f b) { return a + b; }
M3f sub_m3f(M3f a, M3f b) { return a - b; }
M3f mul_m3f(M3f a, M3f b) { return a * b; }
float abs_m3f(M3f a) { return a.norm(); }
V3f mul_m3f_v3f(M3f a, V3f b) { return a * b; }
V3f mul_x3f_v3f(X3f a, V3f b) { return a * b; }

V3f mul_v3f_f(V3f v, float f) { return v * f; }
V3f div_v3f_f(V3f v, float f) { return v / f; }
V3f mul_f_v3f(float f, V3f v) { return f * v; }

M3f mul_m3f_f(M3f m, float f) { return m * f; }
M3f div_m3f_f(M3f m, float f) { return m / f; }
M3f mul_f_m3f(float f, M3f m) { return f * m; }

struct test {
  V3f a;
  float f;
  V3f b;
  int i;
};
namespace std {
template <>
struct is_pod<test> : public std::integral_constant<bool, true> {};
}

void RIFLIB_PYBIND_eigen_types(py::module& m) {
  py::detail::npy_format_descriptor<V3f>::register_dtype();
  py::detail::npy_format_descriptor<M3f>::register_dtype();
  py::detail::npy_format_descriptor<V3<int32_t>>::register_dtype();
  PYBIND11_NUMPY_DTYPE(test, a, i, f, b);

  m.attr("v3f_t") = py::dtype::of<V3f>();
  m.attr("m3f_t") = py::dtype::of<M3f>();
  m.attr("v3i_t") = py::dtype::of<V3<int32_t>>();
  m.attr("test_t") = py::dtype::of<test>();

  m.def("abs_v3f", py::vectorize(abs_v3f));
  m.def("abs_m3f", py::vectorize(abs_m3f));
  m.def("add_v3f", py::vectorize(add_v3f));
  m.def("sub_v3f", py::vectorize(sub_v3f));
  m.def("mul_v3f", py::vectorize(mul_v3f));
  m.def("add_m3f", py::vectorize(add_m3f));
  m.def("sub_m3f", py::vectorize(sub_m3f));
  m.def("mul_m3f", py::vectorize(mul_m3f));
  m.def("mul_m3f_v3f", py::vectorize(mul_m3f_v3f));

  m.def("mul_v3f_f", py::vectorize(mul_v3f_f));
  m.def("div_v3f_f", py::vectorize(div_v3f_f));
  m.def("mul_f_v3f", py::vectorize(mul_f_v3f));
  m.def("mul_m3f_f", py::vectorize(mul_m3f_f));
  m.def("div_m3f_f", py::vectorize(div_m3f_f));
  m.def("mul_f_m3f", py::vectorize(mul_f_m3f));
}
