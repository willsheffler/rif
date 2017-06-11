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
  using R = Ray<F>;
  using T = RayToRay4dHash<R, K>;
  std::string name = "RayToRay4dHash_" + short_str<F>() + short_str<K>();
  py::class_<T>(m, name.c_str())
      .def(py::init<F, F, F>(), "resl"_a = 0.25, "lever"_a = 2, "bound"_a = 128)
      .def_readonly("resl", &T::resl_)
      .def_readonly("lever", &T::lever_)
      .def_readonly("bound", &T::bound_)
      .def_property_readonly("ncells", &T::size)
      .def("get_key", &T::get_key, "ray1"_a, "ray2"_a)
      .def("get_key_aligned", &T::get_key_aligned, "ray"_a)
      .def("get_keys",
           [](T const& self, py::array_t<R> rays1, py::array_t<R> rays2) {
             auto func = [&self](R a, R b) { return self.get_key(a, b); };
             return py::vectorize(func)(rays1, rays2);
           },
           "rays1"_a, "rays2"_a)
      .def("get_center", &T::get_center, "key"_a)
      .def("get_centers",
           [](T const& self, py::array_t<K> keys) {
             auto func = [&self](K key) { return self.get_center(key); };
             return py::vectorize(func)(keys);
           },
           "rays"_a)
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
  using R = Ray<F>;
  using T = Ray5dHash<Ray<F>, K>;
  std::string name = "Ray5dHash_" + short_str<F>() + short_str<K>();
  py::class_<T>(m, name.c_str())
      .def(py::init<F, F, F>(), "resl"_a = 0.25, "lever"_a = 2, "bound"_a = 128)
      .def_readonly("resl", &T::resl_)
      .def_readonly("lever", &T::lever_)
      .def_readonly("bound", &T::bound_)
      .def_property_readonly("ncells", &T::size)
      .def("get_key", &T::get_key, "ray"_a)
      .def("get_keys",
           [](T const& self, py::array_t<R> rays) {
             auto func = [&self](R a) { return self.get_key(a); };
             return py::vectorize(func)(rays);
           },
           "rays"_a)
      .def("get_center", &T::get_center, "key"_a)
      .def("get_centers",
           [](T const& self, py::array_t<K> keys) {
             auto func = [&self](K key) { return self.get_center(key); };
             return py::vectorize(func)(keys);
           },
           "rays"_a)
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
  using R = Ray<F>;
  using T = RayRay10dHash<R, K>;
  std::string name = "RayRay10dHash_" + short_str<F>() + short_str<K>();
  py::class_<T>(m, name.c_str())
      .def(py::init<F, F, F>(), "resl"_a = 0.25, "lever"_a = 2, "bound"_a = 32)
      .def_property_readonly("ncells", &T::size)
      .def("__repr__",
           [=](T const& o) {
             return name + "(resl=" + str(o.resl_) +
                    ", lever=" + str(o.lever_) + ", bound=" + str(o.bound_) +
                    ")";
           })
      .def_readonly("resl", &T::resl_)
      .def_readonly("lever", &T::lever_)
      .def_readonly("bound", &T::bound_)
      .def("get_key", &T::get_key, "ray1"_a, "ray2"_a)
      .def("get_keys",
           [](T const& self, py::array_t<R> rays1, py::array_t<R> rays2) {
             auto func = [&self](R a, R b) { return self.get_key(a, b); };
             return py::vectorize(func)(rays1, rays2);
           },
           "rays1"_a, "rays2"_a)
      .def("get_center", &T::get_center, "key"_a)
      .def("get_centers",
           [](T const& self, py::array_t<K> keys) {
             size_t N = keys.size();
             py::array_t<R> a(N), b(N);
             auto kptr = (K*)keys.request().ptr;
             auto aptr = (R*)a.request().ptr;
             auto bptr = (R*)b.request().ptr;
             for (size_t i = 0; i < N; ++i) {
               auto raypair = self.get_center(kptr[i]);
               aptr[i] = raypair.first;
               bptr[i] = raypair.second;
             }
             return std::make_pair(a, b);
           },
           "keys"_a)

      /* end class_ */;
}

void RIFLIB_PYBIND_geom_ray_hash(py::module& m) {
  m.def("align_ray_pair", &rif::geom::align_ray_pair<float>, "ray1"_a,
        "ray2"_a);
  m.def("align_ray_pairs", py::vectorize(rif::geom::align_ray_pair<float>),
        "ray1"_a, "ray2"_a);
  pybind_ray_hash4<float, int64_t>(m);
  pybind_ray_hash5<float, int64_t>(m);
  pybind_ray_hash10<float, int64_t>(m);
  m.attr("RayToRay4dHash") = m.attr("RayToRay4dHash_f4i8");
  m.attr("Ray5dHash") = m.attr("Ray5dHash_f4i8");
  m.attr("RayRay10dHash") = m.attr("RayRay10dHash_f4i8");
}