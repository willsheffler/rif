#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

#include "rif/actor/Atom.hpp"
#include "rif/eigen_types.hpp"
#include "rif/index/OneSide3dIndex.hpp"

namespace py = pybind11;

template <class XYZ>
void bind_one_sided_index(py::module& m, std::string name) {
  using OSI = rif::index::OneSide3dIndex<XYZ>;
  using F = typename XYZ::Scalar;
  py::class_<OSI>(m, name.c_str())
      .def("__init__",
           [](OSI& instance, py::array_t<XYZ> const& xyzs, F width) {
             new (&instance) OSI(width, xyzs.data(), xyzs.size());
           })
      /**/;
}

void RIFLIB_PYBIND_rif_index(py::module& m) {
  using Atom = rif::actor::Atom<rif::V3f>;
  bind_one_sided_index<Atom>(m, "AtomIndexOneSided");
}
