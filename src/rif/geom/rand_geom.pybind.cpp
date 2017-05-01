#include "pyutil/pybind_numpy.hpp"

#include "geom/rand_geom.hpp"

namespace py = pybind11;
using namespace rif::geom;

void RIFLIB_PYBIND_geom_rand(py::module &m) {
  m.def("rand_normal", py::vectorize(rand_normal<float>));
}
