#pragma once

#include <Eigen/Geometry>

namespace rif {

template <class F>
using V3 = Eigen::Matrix<F, 3, 1>;
template <class F>
using M3 = Eigen::Matrix<F, 3, 3>;
template <class F>
using X3 = Eigen::Transform<F, 3, Eigen::Affine>;

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

using V = Eigen::Vector3f;
using M = Eigen::Matrix3f;
using X = Eigen::Transform<float, 3, Eigen::Affine>;
}
