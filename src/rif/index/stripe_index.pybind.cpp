#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pyutil/pybind_numpy.hpp>

#include "rif/actor/Atom.hpp"
#include "rif/eigen_types.hpp"
#include "rif/index/stripe_index_3d.hpp"

namespace py = pybind11;
using namespace rif;

template <class Index>
void bind_rest(py::class_<Index>& cls) {
  using Pt = typename Index::Pt;
  cls.def("size", &Index::size);
  cls.def("_raw", [](Index const& self, size_t i) { return self.values_[i]; });
  cls.def("nbcount", &Index::template nbcount<Pt>);
}

template <class Point>
void bind_one_sided_index(py::module& m, std::string name) {
  using Index = rif::index::stripe_index_3d<Point>;
  using F = typename Point::Scalar;
  py::class_<Index> cls(m, name.c_str(), py::dynamic_attr());
  cls.def("__init__", [](Index& instance, F width,
                         py::array_t<Point> const& points, py::none) {
    new (&instance) Index();
    instance.init(width, points.data(), points.size());
  });
  bind_rest<Index>(cls);
}

// do not do this, instead if pyobject is being stored,
// keep ref to original payload list / array in Index object,
// store 1->N index into list / array in Index and
// decorate appropriate lookup funcs so they return orig[i];
// also write python factory func to dispatch to appropriate
// Index class and handle the 1->N /
// payload switcharoo decorated functions;

template <class Point, class Payload>
void bind_one_sided_index(py::module& m, std::string name) {
  using Index = rif::index::stripe_index_3d<Point, Payload>;
  using F = typename Point::Scalar;
  py::class_<Index> cls(m, name.c_str(), py::dynamic_attr());
  cls.def("__init__",
          [](Index& instance, F width, py::array_t<Point> const& points,
             py::array_t<Payload> const& payload) {
            if (points.size() != payload.size())
              throw std::invalid_argument("len(points) != len(payload)");
            new (&instance) Index();
            instance.init(width, points.data(), points.size(), payload.data());
          });
  bind_rest<Index>(cls);
}

void RIFLIB_PYBIND_rif_index(py::module& m) {
  using Atom = rif::actor::Atom<V3f>;
  bind_one_sided_index<V3f>(m, "stripe_index_3d_V3");
  bind_one_sided_index<Atom>(m, "stripe_index_3d_Atom");
  bind_one_sided_index<V3f, size_t>(m, "stripe_index_3d_V3_object");
  bind_one_sided_index<Atom, size_t>(m, "stripe_index_3d_Atom_object");
}
