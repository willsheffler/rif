#pragma once

#include <Eigen/Geometry>
#include <iostream>

namespace rif {

template <class F>
using V3 = Eigen::Matrix<F, 3, 1>;
template <class F>
using M3 = Eigen::Matrix<F, 3, 3, Eigen::RowMajor>;  // to match numpy (??)
template <class F>
using X3 = Eigen::Transform<F, 3, Eigen::Affine>;

template <class S>
using A1 = Eigen::Array<S, 1, 1>;
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
template <class S>
using A7 = Eigen::Array<S, 7, 1>;
template <class S>
using A8 = Eigen::Array<S, 8, 1>;
template <class S>
using A9 = Eigen::Array<S, 9, 1>;
template <class S>
using A10 = Eigen::Array<S, 10, 1>;

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

using A1f = A1<float>;
using A2f = A2<float>;
using A3f = A3<float>;
using A4f = A4<float>;
using A5f = A5<float>;
using A6f = A6<float>;
using A7f = A7<float>;
using A8f = A8<float>;
using A9f = A9<float>;
using A10f = A10<float>;

using I1 = A1<int64_t>;
using I2 = A2<int64_t>;
using I3 = A3<int64_t>;
using I4 = A4<int64_t>;
using I5 = A5<int64_t>;
using I6 = A6<int64_t>;
using I7 = A7<int64_t>;
using I8 = A8<int64_t>;
using I9 = A9<int64_t>;
using I10 = A10<int64_t>;

using V3f = V3<float>;
using M3f = M3<float>;
using X3f = X3<float>;
using V3d = V3<double>;
using M3d = M3<double>;
using X3d = X3<double>;

static I1 makeI1(uint64_t a) {
  I1 r;
  r << a;
  return r;
}
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
static I7 makeI7(uint64_t a, uint64_t b, uint64_t c, uint64_t d, uint64_t e,
                 uint64_t f, uint64_t g) {
  I7 r;
  r << a, b, c, d, e, f, g;
  return r;
}
static I8 makeI8(uint64_t a, uint64_t b, uint64_t c, uint64_t d, uint64_t e,
                 uint64_t f, uint64_t g, uint64_t h) {
  I8 r;
  r << a, b, c, d, e, f, g, h;
  return r;
}
static I9 makeI9(uint64_t a, uint64_t b, uint64_t c, uint64_t d, uint64_t e,
                 uint64_t f, uint64_t g, uint64_t h, uint64_t i) {
  I9 r;
  r << a, b, c, d, e, f, g, h, i;
  return r;
}
static I10 makeI10(uint64_t a, uint64_t b, uint64_t c, uint64_t d, uint64_t e,
                   uint64_t f, uint64_t g, uint64_t h, uint64_t i, uint64_t j) {
  I10 r;
  r << a, b, c, d, e, f, g, h, i, j;
  return r;
}

template <class F>
F epsilon2() {
  return std::sqrt(std::numeric_limits<F>::epsilon());
}

template <class F>
std::ostream& operator<<(std::ostream& out, X3<F> x) {
  out << x.rotation() << std::endl;
  out << x.translation().transpose();
  return out;
}
}
