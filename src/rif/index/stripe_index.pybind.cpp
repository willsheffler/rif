#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pyutil/pybind_numpy.hpp>
#include "pyutil/conversions.hpp"

#include "rif/actor/Atom.hpp"
#include "rif/eigen_types.hpp"
#include "rif/index/stripe_index_3d.hpp"

namespace py = pybind11;
using namespace rif;
using namespace rif::util;
using namespace pyutil;

template <class Point>
bool is_simple_point() {
  return std::is_same<Point, V3f>::value || std::is_same<Point, V3d>::value;
}

template <class Index>
void bind_rest(py::class_<Index>& cls) {
  using Point = typename Index::Pt;
  using F = typename Point::Scalar;
  cls.def("__len__", &Index::size);
  cls.def("_raw", [](Index const& self, size_t i) { return self.values_[i]; });
  cls.def("_raw_point", [](Index const& self, size_t i) {
    return get_first_if_pair(self.values_[i]);
  });
  cls.def("neighbor_count", &Index::template nbcount<Point>);
  // if (is_simple_point<Point>()) {
  // cls.def("neighbor_count", [](Index const& self, py::object obj) {
  // return self.nbcount(sequence_to<Point, F>(obj, 3));
  // });
  // }
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
  cls.def("neighbors", [](Index const& self, Point query) {
    return pyutil::to_pyarray(self.neighboring_points(query));
  });
  if (is_simple_point<Point>())
    cls.def("neighbors", [](Index const& self, py::object obj) {
      Point query = pyutil::sequence_to<Point, F>(obj, 3);
      return pyutil::to_pyarray(self.neighboring_points(query));
    });
  bind_rest<Index>(cls);
}

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
  cls.def("_raw_payload", [](Index const& self, size_t i) {
    return get_second_if_pair(self.values_[i]);
  });
  cls.def("neighbors", [](Index const& self, Point query) {
    return pyutil::to_pyarray(self.neighboring_payloads(query));
  });
  if (is_simple_point<Point>())
    cls.def("neighbors", [](Index const& self, py::object obj) {
      Point query = pyutil::sequence_to<Point, float>(obj, 3);
      return pyutil::to_pyarray(self.neighboring_payloads(query));
    });
  bind_rest<Index>(cls);
}

void RIFLIB_PYBIND_rif_index(py::module& m) {
  using Atom = rif::actor::Atom<V3f>;
  bind_one_sided_index<V3f>(m, "stripe_index_3d_V3");
  bind_one_sided_index<Atom>(m, "stripe_index_3d_Atom");
  bind_one_sided_index<V3f, V3f>(m, "stripe_index_3d_V3_V3");
  bind_one_sided_index<V3f, Atom>(m, "stripe_index_3d_V3_Atom");
  bind_one_sided_index<V3f, size_t>(m, "stripe_index_3d_V3_object");
  bind_one_sided_index<Atom, size_t>(m, "stripe_index_3d_Atom_object");
  m.def("make_huge_array", [](size_t s) {
    int* raw = (int*)malloc(s * sizeof(int));
    raw[0] = 293847;
    std::cout << "done alloc size " << s << std::endl;
    for (int i = 0; i < 10; ++i) std::cout << raw[i] << std::endl;
    auto a = py::array(pybind11::buffer_info(
        raw, sizeof(int), pybind11::format_descriptor<int>::format(), 1, {s},
        {sizeof(int)}));
    std::cout << "done making py::array" << std::endl;
    return a;

  });
}
