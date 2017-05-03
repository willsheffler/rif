#include "pyutil/pybind_numpy.hpp"

#include <rif/geom/Ray.hpp>

namespace py = pybind11;
using namespace rif::geom;
using Rf = Ray<float>;
using Rd = Ray<double>;

namespace std {
template <class F>
struct is_pod<Ray<F>> : public std::integral_constant<bool, true> {};
}

void RIFLIB_PYBIND_geom_Ray(py::module &m) {
  PYBIND11_NUMPY_DTYPE(Rf, orig, dirn);
  PYBIND11_NUMPY_DTYPE(Rd, orig, dirn);
  m.def("rand_ray_gaussian", py::vectorize(rand_ray_gaussian<float>));
  py::class_<Rd>(m, "Ray")
      .def_readwrite("orig", &Rd::orig)
      .def_readwrite("dirn", &Rd::dirn)
      .def("__getitem__",
           [](Rd &self, int i) {
             if (0 == i) return self.orig;
             if (1 == i) return self.dirn;
             throw std::invalid_argument("Index out of range");
             return self.orig;
           })
      .def("__repr__",
           [](Rd const &rd) {
             std::ostringstream oss;
             oss << "Ray(orig=[" << rd.orig[0] << ", " << rd.orig[1] << ", "
                 << rd.orig[2] << "], dirn=[" << rd.dirn[0] << ", "
                 << rd.dirn[1] << ", " << rd.dirn[2] << "])";
             return oss.str();
           })
      /**/;
}
