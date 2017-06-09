#include "rif/geom/ray_hash.hpp"
#include "rif/util/str.hpp"

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

namespace py = pybind11;
using namespace pybind11::literals;
using namespace rif::geom;
using rif::util::short_str;
using rif::util::str;

template <class F, class K>
void pybind_ray_hash4(py::module& m) {
  using T = RayToRay4dHash<Ray<F>, K>;
  std::string name = "RayToRay4dHash_" + short_str<F>() + short_str<K>();
  py::class_<T>(m, name.c_str())
      .def(py::init<F, F, F>(), "resl"_a = 0.25, "lever"_a = 2, "bound"_a = 32)
      .def("__repr__",
           [=](T const& o) {
             return name + "(resl=" + str(o.resl_) +
                    ", lever=" + str(o.lever_) + ", bound=" + str(o.bound_) +
                    ")";
           })
      /* end class_ */;
}

template <class F, class K>
void pybind_ray_hash5(py::module& m) {
  using T = Ray5dHash<Ray<F>, K>;
  std::string name = "Ray5dHash_" + short_str<F>() + short_str<K>();
  py::class_<T>(m, name.c_str())
      .def(py::init<F, F, F>(), "resl"_a = 0.25, "lever"_a = 2, "bound"_a = 32)
      .def("__repr__",
           [=](T const& o) {
             return name + "(resl=" + str(o.resl_) +
                    ", lever=" + str(o.lever_) + ", bound=" + str(o.bound_) +
                    ")";
           })
      /* end class_ */;
}

template <class F, class K>
void pybind_ray_hash10(py::module& m) {
  using T = RayRay10dHash<Ray<F>, K>;
  std::string name = "RayRay10dHash_" + short_str<F>() + short_str<K>();
  py::class_<T>(m, name.c_str())
      .def(py::init<F, F, F>(), "resl"_a = 0.25, "lever"_a = 2, "bound"_a = 32)
      .def("__repr__",
           [=](T const& o) {
             return name + "(resl=" + str(o.resl_) +
                    ", lever=" + str(o.lever_) + ", bound=" + str(o.bound_) +
                    ")";
           })
      .def_readonly("resl", &T::resl_)
      .def_readonly("lever", &T::lever_)
      .def_readonly("bound", &T::bound_)
      // .def("get_key", py::vectorize<T::get_key>)
      .def("get_key", &T::get_key, "ray1"_a, "ray2"_a)
      .def("get_center", &T::get_center, "key"_a)
      /* end class_ */;
}

void RIFLIB_PYBIND_geom_ray_hash(py::module& m) {
  pybind_ray_hash4<float, int64_t>(m);
  pybind_ray_hash5<float, int64_t>(m);
  pybind_ray_hash10<float, int64_t>(m);
  m.attr("RayToRay4dHash") = m.attr("RayToRay4dHash_f4i8");
  m.attr("Ray5dHash") = m.attr("Ray5dHash_f4i8");
  m.attr("RayRay10dHash") = m.attr("RayRay10dHash_f4i8");
}