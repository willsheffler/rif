#pragma once

#include <Eigen/Geometry>

namespace rif {

template <class F>
using V3 = Eigen::Matrix<F, 3, 1>;
template <class F>
using M3 = Eigen::Matrix<F, 3, 3, Eigen::RowMajor>;  // to match numpy (??)
template <class F>
using X3 = Eigen::Transform<F, 3, Eigen::Affine>;

template <class S>
using A2 = Eigen::Array<S, 2, 1>;
template <class S>
using A3 = Eigen::Array<S, 3, 1>;
template <class S>
using A4 = Eigen::Array<S, 4, 1>;
template <class S>
using A5 = Eigen::Array<S, 5, 1>;
template <class S>
using A6 = Eigen::Array<S, 6, 1>;

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
using I5 = A5<int64_t>;
using I6 = A6<int64_t>;
using V3f = V3<float>;
using M3f = M3<float>;
using X3f = X3<float>;
using V3d = V3<double>;
using M3d = M3<double>;
using X3d = X3<double>;

static I5 makeI5(uint64_t a, uint64_t b, uint64_t c, uint64_t d, uint64_t e) {
  I5 r;
  r << a, b, c, d, e;
  return r;
}
static I6 makeI6(uint64_t a, uint64_t b, uint64_t c, uint64_t d, uint64_t e,
                 uint64_t f) {
  I6 r;
  r << a, b, c, d, e, f;
  return r;
}
}
