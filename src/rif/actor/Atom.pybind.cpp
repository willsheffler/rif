#include "pybind11/numpy.h"
#include "pybind11/pybind11.h"

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

void RIFLIB_PYBIND_actor_Atom(py::module &m) {
  ::pybind11::detail::npy_format_descriptor<Atom<V3f>>::register_dtype({
      ::pybind11::detail::field_descriptor{
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
          ::pybind11::detail::npy_format_descriptor<float>::dtype()},
      ::pybind11::detail::field_descriptor{
          "atype", 12, sizeof(int16_t), alignof(int16_t),
          ::pybind11::format_descriptor<int16_t>::format(),
          ::pybind11::detail::npy_format_descriptor<int16_t>::dtype()},
      ::pybind11::detail::field_descriptor{
          "rtype", 14, sizeof(int8_t), alignof(int8_t),
          ::pybind11::format_descriptor<int8_t>::format(),
          ::pybind11::detail::npy_format_descriptor<int8_t>::dtype()},
      ::pybind11::detail::field_descriptor{
          "anum", 15, sizeof(int8_t), alignof(int8_t),
          ::pybind11::format_descriptor<int8_t>::format(),
          ::pybind11::detail::npy_format_descriptor<int8_t>::dtype()},
  });

  m.attr("atom_t") = py::dtype::of<Atom<V3f>>();
}
