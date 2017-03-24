#ifndef INCLUDED_numeric_ray_ray_hash_HH
#define INCLUDED_numeric_ray_ray_hash_HH

#include <iostream>
#include <numeric/Ray.hpp>
#include <numeric/types.hpp>

namespace rif {
namespace numeric {

///@brief
///@detail

/**
 * @brief       Align rays to ray a's canonical position
 *
 * @param[in]  a     ray1
 * @param[in]  b     ray2
 *
 * @tparam     F     float type
 *
 * @return     Ray b aligned to Ray a's canonical position
 * @detail     canonical position places r at the origin along x and s in the
 * x,y plane
 */
template <class F>
Ray<F> align_ray_pair(Ray<F> a, Ray<F> b) {
  M3<F> rotation;
  V3<F> basis1 = a.direction;
  V3<F> basis3 = basis1.cross(b.origin - a.origin).normalized();
  V3<F> basis2 = basis3.cross(basis1).normalized();
  rotation.row(0) = basis1;
  rotation.row(1) = basis2;
  rotation.row(2) = basis3;
  V3<F> a2_origin = rotation * a.origin;
  Ray<F> b2 = rotation * b;
  b2.origin -= a2_origin;
  // std::cout << rotation.determinant() << std::endl;
  // Ray<F> a2 = rotation * a;
  // std::cout << "a  " << a << std::endl;
  // std::cout << "b  " << b << std::endl;
  // std::cout << "a2 " << a2 << std::endl;
  // std::cout << "b2 " << b2 << std::endl;
  return b2;
}
}
}

#endif
