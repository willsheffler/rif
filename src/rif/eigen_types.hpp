#pragma once

#include <Eigen/Geometry>
#include <iostream>

namespace rif {

template <class F>
using V3 = Eigen::Matrix<F, 3, 1>;  // todo: should use aligned vector3?
template <class F>
using M3 = Eigen::Matrix<F, 3, 3, Eigen::RowMajor>;  // to match numpy (??)
template <class F>
using X3 = Eigen::Transform<F, 3, Eigen::Affine, Eigen::RowMajor>;
template <class F>
using Xc3 = Eigen::Transform<F, 3, Eigen::AffineCompact>;  // always colmajor

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
using Xc3f = Xc3<float>;
using V3d = V3<double>;
using M3d = M3<double>;
using X3d = X3<double>;
using Xc3d = Xc3<double>;

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

static_assert(sizeof(V3f) == 4 * 3, "V3f insane size!");
static_assert(sizeof(M3f) == 4 * 3 * 3, "M3f insane size!");
static_assert(sizeof(X3f) == 4 * 4 * 4, "X3f insane size");
static_assert(sizeof(Xc3f) == 4 * 3 * 4, "X3f insane size");

template <class X>
using ScalarOf = typename X::Scalar;

namespace util {
using namespace Eigen;

// this shit is dumb... use ::linear and ::translation instead

template <class F>
size_t raw_asif_rowmajor(Transform<F, 3, Affine, RowMajor>& x, size_t i) {
  static size_t idx[12] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};
  return x.data()[idx[i]];
}

template <class F>
size_t raw_asif_rowmajor(Transform<F, 3, Affine, ColMajor>& x, size_t i) {
  static size_t idx[12] = {0, 4, 8, 12, 1, 5, 9, 13, 2, 6, 10, 14};
  return x.data()[idx[i]];
}

// AffineCompact is always ColMajor!?!
template <class F>
size_t raw_asif_rowmajor(Transform<F, 3, AffineCompact>& x, size_t i) {
  static size_t idx[12] = {0, 3, 6, 9, 1, 4, 7, 10, 2, 5, 8, 11};
  return x.data()[idx[i]];
}
}  // namespace util

}  // namespace rif

namespace std {
template <class F>
struct is_pod<rif::V3<F>> : public std::integral_constant<bool, true> {};
template <class F>
struct is_pod<rif::M3<F>> : public std::integral_constant<bool, true> {};
template <class F>
struct is_pod<rif::X3<F>> : public std::integral_constant<bool, true> {};
}