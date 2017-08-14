#include <pybind11/numpy.h>
#include <pybind11/operators.h>
#include "pyutil/pybind_numpy.hpp"

#include "eigen_types.hpp"

#include <iostream>

namespace py = pybind11;
using namespace rif;

template <class A, class B, class C>
C add(A a, B b) {
  return a + b;
}
template <class A, class B, class C>
C sub(A a, B b) {
  return a - b;
}
template <class A, class B, class C>
C mul(A a, B b) {
  return a * b;
}
template <class A, class B, class C>
C div(A a, B b) {
  return a / b;
}

V3f add_v3f(V3f a, V3f b) { return a + b; }
V3f sub_v3f(V3f a, V3f b) { return a - b; }
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

template <class M>
void bind_eigen_matrix_fixed(py::module& m, std::string name) {
  using F = ScalarOf<M>;
  auto cls = py::class_<M>(m, name.c_str(), py::buffer_protocol());
  cls.def("__eq__", [](M const& m, M const& n) { return m.isApprox(n); })
      .def("__getitem__", [](M const& m, int i) { return m.data()[i]; })
      .def("__setitem__", [](M& m, int i, F f) { return m.data()[i] = f; })
      .def("__len__", [](M const& m) { return m.rows() * m.cols(); })
      .def("__abs__", [](M const& m) { return m.norm(); })
      .def("__add__", [](M const& a, M const& b) { return M(a + b); },
           py::is_operator())
      .def("__sub__", [](M const& a, M const& b) { return M(a - b); },
           py::is_operator())
      .def_property_readonly("norm", [](M const& m) { return m.norm(); })
      .def("normalize", [](M& m) { m.normalize(); })
      .def_property_readonly("dtype",
                             [](M const&) { return py::dtype::of<M>(); })
      .def_property_readonly_static(
          "dtype", [](py::object) { return py::dtype::of<M>(); })
      .def_buffer([](M& m) -> py::buffer_info {
        std::vector<size_t> shape = {(size_t)m.rows(), (size_t)m.cols()};
        std::vector<size_t> stride = {sizeof(F) * m.rows(), sizeof(F)};
        if (m.cols() == 1) {
          shape.pop_back();
          stride.pop_back();
          stride[0] = sizeof(F);
        }
        return py::buffer_info(m.data(), sizeof(F),
                               py::format_descriptor<F>::format(), shape.size(),
                               shape, stride);
      })
      /**/;
  py::detail::npy_format_descriptor<M>::register_dtype();
  // TODO: somehow set dtype.type here?!?
  // quaternion_descr->typeobj = &PyQuaternion_Type;
}

template <class X>
void bind_eigen_xform(py::module& m, std::string name) {
  auto cls = py::class_<X>(m, name.c_str() /*, buffer?*/);
  cls.def_property_readonly_static(
      "dtype", [](py::object) { return py::dtype::of<X>(); })
      /**/;
  py::detail::npy_format_descriptor<X>::register_dtype();
}

void print_numpy_info(py::array a) {
  py::buffer_info b(py::buffer(a).request());
  std::cout << "size: " << b.size << std::endl;
  std::cout << "itemsize: " << b.itemsize << std::endl;
  std::cout << "ndim: " << b.ndim << std::endl;
  std::cout << "format: '" << b.format << "'" << std::endl;
}

void RIFLIB_PYBIND_eigen_types(py::module& m) {
  bind_eigen_matrix_fixed<V3f>(m, "V3");
  bind_eigen_matrix_fixed<V3d>(m, "V3d");
  bind_eigen_matrix_fixed<V3<int32_t>>(m, "V3i4");
  bind_eigen_matrix_fixed<V3<int64_t>>(m, "V3i");
  bind_eigen_matrix_fixed<M3f>(m, "M3");
  bind_eigen_matrix_fixed<M3d>(m, "M3d");
  bind_eigen_xform<X3f>(m, "X3");
  bind_eigen_xform<X3d>(m, "X3d");
  bind_eigen_xform<Xc3f>(m, "Xc3");
  bind_eigen_xform<Xc3d>(m, "Xc3d");

  /*
    for o in "add sub mul div".split():
        for t1 in "float V3f M3f X3f".split():
            for t2 in "float V3f M3f X3f".split():
                r = 'RET'
                print(('m.def("%(o)s_%(t1)s_%(t2)s", ' +
                       'py::vectorize(%(o)s<%(t1)s, %(t2)s, %(r)s>));') %
    vars())

  */
  m.def("rifop_abs_V3", py::vectorize(abs_v3f));
  m.def("rifop_abs_M3", py::vectorize(abs_m3f));
  m.def("rifop_add_V3_V3", py::vectorize(add<V3f, V3f, V3f>));
  m.def("rifop_add_M3_M3", py::vectorize(add<M3f, M3f, M3f>));
  m.def("rifop_sub_V3_V3", py::vectorize(sub<V3f, V3f, V3f>));
  m.def("rifop_sub_M3_M3", py::vectorize(sub<M3f, M3f, M3f>));
  m.def("rifop_mul_fl_V3", py::vectorize(mul<float, V3f, V3f>));
  m.def("rifop_mul_fl_M3", py::vectorize(mul<float, M3f, M3f>));
  m.def("rifop_mul_V3_fl", py::vectorize(mul<V3f, float, V3f>));
  m.def("rifop_mul_M3_V3", py::vectorize(mul<M3f, V3f, V3f>));
  m.def("rifop_mul_M3_M3", py::vectorize(mul<M3f, M3f, M3f>));
  m.def("rifop_mul_X3_V3", py::vectorize(mul<X3f, V3f, V3f>));
  m.def("rifop_mul_X3_M3", py::vectorize(mul<X3f, M3f, M3f>));
  m.def("rifop_mul_X3_X3", py::vectorize(mul<X3f, X3f, X3f>));
  m.def("rifop_div_V3_fl", py::vectorize(div<V3f, float, V3f>));
  m.def("rifop_div_M3_fl", py::vectorize(div<M3f, float, M3f>));

  // testing crap
  m.def("print_numpy_info", &print_numpy_info);
  // PYBIND11_NUMPY_DTYPE(test, a, i, f, b);
  // m.attr("test_t") = py::dtype::of<test>();
}
