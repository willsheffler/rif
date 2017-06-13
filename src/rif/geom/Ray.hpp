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
  using VM = Eigen::Map<V3<F>, 0, Eigen::InnerStride<2>>;
  using VMconst = Eigen::Map<V3<F> const, 0, Eigen::InnerStride<2>>;
  using M42 = Eigen::Matrix<F, 4, 2, Eigen::RowMajor>;
  M42 _m42;
  Ray() { _m42 << 0, 1, 0, 0, 0, 0, 1, 0; }
  Ray(V3<F> o, V3<F> d) {
    orig() = o;
    dirn() = d.normalized();
    _m42.row(3) = A2<F>(1, 0);  // homo-coords, 1==point, 0=dirn
  }
  template <class F2, int OPTS>
  Ray(Eigen::Matrix<F2, 4, 2, OPTS> const& in) : _m42(in) {}
  VM orig() { return VM(_m42.data()); }
  VM dirn() { return VM(_m42.data() + 1); }
  VMconst orig() const { return VMconst(_m42.data()); }
  VMconst dirn() const { return VMconst(_m42.data() + 1); }
  void normalize() { dirn().normalize(); }
  V getorig() { return orig(); }
  V getdirn() { return dirn(); }
  void setorig(V o) { orig() = o; }
  void setdirn(V d) { dirn() = d.normalized(); }
};

template <class F>
std::ostream& operator<<(std::ostream& s, Ray<F> r) {
  s << "Ray( " << r.orig().transpose() << " ; " << r.dirn().transpose() << " )";
  return s;
}

template <class F1, class F2>
Ray<F2> operator*(X3<F1> const& x, Ray<F2> const& r) {
  using namespace Eigen;
  return Ray<F2>(typename Ray<F2>::M42(x.matrix() * r._m42));

  // return Ray<F2>(x * r.orig(), x.linear() * r.dirn());
}
template <class F1, class F2>
Ray<F2> operator*(M3<F1> const& m, Ray<F2> const& r) {
  return Ray<F2>(m * r.orig(), m * r.dirn());
  // return Ray<F2>(m * r._m42);
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
  r.dirn()[0] = gaussian(rng);
  r.dirn()[1] = gaussian(rng);
  r.dirn()[2] = gaussian(rng);
  r.dirn().normalize();
  r.orig()[0] = gaussian(rng) * sd;
  r.orig()[1] = gaussian(rng) * sd;
  r.orig()[2] = gaussian(rng) * sd;
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
