#include "pybind11/numpy.h"
#include "pybind11/pybind11.h"

#include "eigen_types.hpp"

namespace py = pybind11;

namespace std {
template <class F, int M, int N>
struct is_pod<Eigen::Matrix<F, M, N>>
    : public std::integral_constant<bool, true> {};
}

struct fxyz {
  float x;
  float y;
  float z;
};

void RIFLIB_PYBIND_eigen_types(py::module &m) {
  PYBIND11_NUMPY_DTYPE(fxyz, x, y, z);
  m.attr("fxyz_t") = py::dtype::of<fxyz>();
}
