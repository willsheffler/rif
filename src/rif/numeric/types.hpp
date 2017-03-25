#pragma once

#include <Eigen/Geometry>

namespace rif {

template <class F>
using V3 = Eigen::Matrix<F, 3, 1>;
template <class F>
using M3 = Eigen::Matrix<F, 3, 3>;
template <class F>
using X3 = Eigen::Transform<F, 3, Eigen::Affine>;

template <class S>
using A2 = Eigen::Array<S, 2, 1>;
template <class S>
using A3 = Eigen::Array<S, 3, 1>;
template <class S>
using A4 = Eigen::Array<S, 4, 1>;

template <class F>
X3<F> xform(M3<F> m, V3<F> v) {
  X3<F> x(m);
  x.translation() = v;
  return x;
}
template <class F>
X3<F> xform(V3<F> v) {
  M3<F> m = M3<F>::Identity();
  return xform(m, v);
}

using I2 = A2<int64_t>;
using I3 = A3<int64_t>;
using I4 = A4<int64_t>;
using V3f = V3<float>;
using M3f = M3<float>;
using X3f = X3<float>;
using V3d = V3<double>;
using M3d = M3<double>;
using X3d = X3<double>;
}
