#pragma once

#include <Eigen/Dense>
#include <iostream>
#include <numeric/types.hpp>
#include <random>

namespace rif {
namespace numeric {

/**
 * @brief      a ray in space
 *
 * @tparam     F     float type
 */
template <class F = float>
struct Ray {
  using Scalar = F;
  using V = V3<F>;
  V3<F> origin;
  V3<F> direction;
  Ray() : origin(0, 0, 0), direction(1, 0, 0) {}
  Ray(V3<F> o, V3<F> d) : origin(o), direction(d.normalized()) {}
};

template <class F>
std::ostream& operator<<(std::ostream& s, Ray<F> r) {
  s << "Ray( " << r.origin.transpose() << " ; " << r.direction.transpose()
    << " )";
  return s;
}

template <class F1, class F2>
Ray<F2> operator*(X3<F1> const& x, Ray<F2> const& r) {
  return Ray<F2>(x * r.origin, x.rotation() * r.direction);
}
template <class F1, class F2>
Ray<F2> operator*(M3<F1> const& m, Ray<F2> const& r) {
  return Ray<F2>(m * r.origin, m * r.direction);
}

/**
 * @brief      generate random ray
 *
 * @param      rng
 * @param[in]  sd    std dev for origin position
 *
 * @tparam     RNG
 * @tparam     F     float type
 *
 * @return     { description_of_the_return_value }
 */
template <class RNG, class F = float>
Ray<F> rand_ray_gaussian(RNG& rng, float sd = 10.0) {
  Ray<F> r;
  std::normal_distribution<> gaussian;
  r.direction[0] = gaussian(rng);
  r.direction[1] = gaussian(rng);
  r.direction[2] = gaussian(rng);
  r.direction.normalize();
  r.origin[0] = gaussian(rng) * sd;
  r.origin[1] = gaussian(rng) * sd;
  r.origin[2] = gaussian(rng) * sd;
  return r;
}
}
}
