#ifndef INCLUDED_numeric_util_HH
#define INCLUDED_numeric_util_HH

#include <math.h>
#include <limits>

namespace rif {
namespace numeric {

template <class X>
using scalar = typename X::Scalar;

template <class Position>
bool approx_eq(Position const &a, Position const &b) {
  return a.isApprox(b);
}

///@brief return sum of highest two elements in vector
template <class Vector, class Index>
void max2(Vector vector, typename Vector::Scalar &mx1,
          typename Vector::Scalar &mx2, Index &argmax_1, Index &argmax_2) {
  // TODO: is there a faster way to do this?
  mx1 = vector.maxCoeff(&argmax_1);
  vector[argmax_1] = -std::numeric_limits<typename Vector::Scalar>::max();
  mx2 = vector.maxCoeff(&argmax_2);
}

template <class Float>
Float rad2quat(Float rad) {
  return sqrt(1.0 - cos(rad / 2.0) * cos(rad / 2.0));
}
template <class Float>
Float deg2quat(Float deg) {
  Float rad = deg * M_PI / 180.0;
  return rad2quat(rad);
}

template <class Float>
Float sqr(Float const &r) {
  return r * r;
}

template <class Float>
Float sigmoidish(Float const &sqdist, Float const &mindis,
                 Float const &maxdis) {
  if (sqdist > maxdis * maxdis) {
    return 0.0;
  } else if (sqdist < mindis * mindis) {
    return 1.0;
  } else {
    Float dist = sqrt(sqdist);
    return sqr(1.0 - sqr((dist - mindis) / (maxdis - mindis)));
  }
}
}
}

#endif
