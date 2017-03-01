#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <sampling/orientations.hpp>

#include <Eigen/Geometry>

namespace py = pybind11;

void RIFLIB_PYBIND_sampling_orientations(py::module &m) {
  m.def("read_karney_orientation_file", &read_karney_orientation_file);
}
