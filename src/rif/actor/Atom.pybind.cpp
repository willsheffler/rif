#include "pyutil/pybind_numpy.hpp"

#include <Eigen/Dense>
#include "pyutil/nonpod_dtype_support.h"
#include "rif/actor/Atom.hpp"
#include "rif/eigen_types.hpp"

namespace py = pybind11;
using namespace rif;
using namespace rif::actor;

using F = float;
using A = Atom<V3f>;
using V = V3f;
using M = M3f;
using X = X3f;

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

float abs_Atom(A a) { return a.pos.norm(); }

void RIFLIB_PYBIND_actor_Atom(py::module &m) {
  static_assert(sizeof(Atom<V3f>) == 16, "bad atom size");
  PYBIND11_NUMPY_DTYPE(A, pos, atype, anum, rtype);
  py::class_<Atom<V3f>> acls(m, "Atom");
  acls.def("__repr__", [](Atom<V3f> const &self) -> std::string {
    std::ostringstream oss;
    oss << self;
    return oss.str();
  });
  acls.attr("dtype") = py::dtype::of<Atom<V3f>>();

  m.def("rifop_abs_AT", py::vectorize(abs_Atom));

  m.def("rifop_add_AT_V3", py::vectorize(add<A, V, A>));
  m.def("rifop_add_V3_AT", py::vectorize(add<V, A, A>));
  m.def("rifop_sub_AT_V3", py::vectorize(sub<A, V, A>));
  m.def("rifop_sub_V3_AT", py::vectorize(sub<V, A, A>));
  m.def("rifop_mul_fl_AT", py::vectorize(mul<F, A, A>));
  m.def("rifop_mul_AT_fl", py::vectorize(mul<A, F, A>));
  m.def("rifop_mul_M3_AT", py::vectorize(mul<M, A, A>));
  m.def("rifop_mul_X3_AT", py::vectorize(mul<X, A, A>));
  m.def("rifop_div_AT_fl", py::vectorize(div<A, F, A>));
}
