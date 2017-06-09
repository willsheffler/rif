#pragma once

#include <eigen_types.hpp>
#include <iostream>
#include "rif/global_rng.hpp"

namespace rif {
namespace geom {

/**
 * @brief      a ray in space
 *
 * @tparam     F     float type
 */
template <class F = float>
struct Ray {
  using Scalar = F;
  using V = V3<F>;
  V3<F> orig;
  V3<F> dirn;
  Ray() : orig(0, 0, 0), dirn(1, 0, 0) {}
  Ray(V3<F> o, V3<F> d) : orig(o), dirn(d.normalized()) {}
};

template <class F>
std::ostream& operator<<(std::ostream& s, Ray<F> r) {
  s << "Ray( " << r.orig.transpose() << " ; " << r.dirn.transpose() << " )";
  return s;
}

template <class F1, class F2>
Ray<F2> operator*(X3<F1> const& x, Ray<F2> const& r) {
  return Ray<F2>(x * r.orig, x.rotation() * r.dirn);
}
template <class F1, class F2>
Ray<F2> operator*(M3<F1> const& m, Ray<F2> const& r) {
  return Ray<F2>(m * r.orig, m * r.dirn);
}

/**
 * @brief      generate random ray
 *
 * @param      rng
 * @param[in]  sd    std dev for orig position
 *
 * @tparam     RNG
 * @tparam     F     float type
 *
 * @return     { description_of_the_return_value }
 */
template <class F = float>
Ray<F> rand_ray_gaussian_rng(std::mt19937& rng, float sd = 10.0) {
  Ray<F> r;
  std::normal_distribution<> gaussian;
  r.dirn[0] = gaussian(rng);
  r.dirn[1] = gaussian(rng);
  r.dirn[2] = gaussian(rng);
  r.dirn.normalize();
  r.orig[0] = gaussian(rng) * sd;
  r.orig[1] = gaussian(rng) * sd;
  r.orig[2] = gaussian(rng) * sd;
  return r;
}
template <class F = float>
Ray<F> rand_ray_gaussian(F sd = 10.0) {
  auto& rng(global_rng());
  return rand_ray_gaussian_rng<F>(rng, sd);
}
}
}

// for pybind to know Ray 'is_pod_struct'
namespace std {
template <class F>
struct is_pod<rif::geom::Ray<F>> : public std::integral_constant<bool, true> {};
template <class F>
struct is_standard_layout<rif::geom::Ray<F>>
    : public std::integral_constant<bool, true> {};
#if !defined(__GNUG__) || defined(__clang__) || __GNUC__ >= 5
template <class F>
struct is_trivially_copyable<rif::geom::Ray<F>>
    : public std::integral_constant<bool, true> {};
#else
// GCC 4 doesn't implement is_trivially_copyable, so approximate it
template <class F>
struct is_trivially_destructible<rif::geom::Ray<F>>
    : public std::integral_constant<bool, true> {};
template <class F>
struct has_trivial_copy_constructor<rif::geom::Ray<F>>
    : public std::integral_constant<bool, true> {};
template <class F>
struct has_trivial_copy_assign<rif::geom::Ray<F>>
    : public std::integral_constant<bool, true> {};
#endif
}
