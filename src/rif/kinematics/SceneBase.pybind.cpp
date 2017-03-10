#include "kinematics/SceneBase.hpp"
#include <pybind11/pybind11.h>
#include <Eigen/Geometry>
namespace py = pybind11;

// 2

struct teststruct {};

void RIFLIB_PYBIND_kinematics_SceneBase(py::module &m) {
  using Xform = Eigen::Transform<float, 3, Eigen::AffineCompact>;
  // py::class_<rif::kinematics::SceneBase<Xform, uint64_t>>(m, "Scene")
  py::class_<teststruct>(m, "Scene").def(py::init<>());
}