#include "pyutil/pybind_numpy.hpp"

#include <rif/geom/Ray.hpp>

namespace py = pybind11;
using namespace rif::geom;
using R = Ray<float>;

namespace std {
template <class F>
struct is_pod<Ray<F>> : public std::integral_constant<bool, true> {};
}

void RIFLIB_PYBIND_geom_Ray(py::module &m) {
  PYBIND11_NUMPY_DTYPE(R, orig, dirn);
  m.def("rand_ray_gaussian", py::vectorize(rand_ray_gaussian<float>));
}
