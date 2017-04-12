#include "pyutil/pybind_numpy.hpp"

#include <Eigen/Dense>
#include "actor/Atom.hpp"
#include "eigen_types.hpp"

namespace py = pybind11;
using namespace rif;
using namespace rif::actor;

namespace std {
template <class P>
struct is_pod<rif::actor::Atom<P>> : public std::integral_constant<bool, true> {
};
}

using A = Atom<V3f>;
using V = V3f;
using M = M3f;
using X = X3f;

A add_v3f_atom(V v, A a) { return A(v + a.pos, a); }
A add_atom_v3f(A a, V v) { return A(a.pos + v, a); }
A sub_v3f_atom(V v, A a) { return A(v - a.pos, a); }
A sub_atom_v3f(A a, V v) { return A(a.pos - v, a); }

void RIFLIB_PYBIND_actor_Atom(py::module &m) {
  static_assert(sizeof(Atom<V3f>) == 16, "bad atom size");
  PYBIND11_NUMPY_DTYPE(A, pos, atype, rtype, anum);
  m.attr("atom_t") = py::dtype::of<Atom<V3f>>();
  m.def("add_v3f_atom", py::vectorize(add_v3f_atom));
  m.def("add_atom_v3f", py::vectorize(add_atom_v3f));
  m.def("sub_v3f_atom", py::vectorize(sub_v3f_atom));
  m.def("sub_atom_v3f", py::vectorize(sub_atom_v3f));
}
