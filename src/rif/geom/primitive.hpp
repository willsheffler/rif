#pragma once

#include <eigen_types.hpp>
#include <iostream>
#include <vector>

namespace rif {
namespace geom {

template <class Ary>
auto welzl_bounding_sphere(Ary const& points) noexcept;

template <class F>
class Sphere {
  using Scalar = F;
  using Vec3 = V3<F>;
  using Mat3 = M3<F>;

 public:
  Vec3 center;
  F radius;

  constexpr static const F epsilon =
      std::sqrt(std::numeric_limits<F>::epsilon());

  Sphere() : center(0, 0, 0), radius(1) {}
  Sphere(Vec3 c, F r) : center(c), radius(r) {}
  Sphere(Vec3 O) : center(O), radius(epsilon) {}
  Sphere(Vec3 O, Vec3 A) {
    Vec3 a = A - O;
    Vec3 o = 0.5 * a;
    radius = o.norm() + epsilon;
    center = O + o;
  }
  Sphere(Vec3 O, Vec3 A, Vec3 B) {
    Vec3 a = A - O, b = B - O;
    F det_2 = 2.0 * ((a.cross(b)).dot(a.cross(b)));
    Vec3 o = (b.dot(b) * ((a.cross(b)).cross(a)) +
              a.dot(a) * (b.cross(a.cross(b)))) /
             det_2;
    radius = o.norm() + epsilon;
    center = O + o;
  }
  Sphere(Vec3 O, Vec3 A, Vec3 B, Vec3 C) {
    Vec3 a = A - O, b = B - O, c = C - O;
    Mat3 cols;
    cols.col(0) = a;
    cols.col(1) = b;
    cols.col(2) = c;
    F det_2 = 2.0 * Mat3(cols).determinant();
    Vec3 o = (c.dot(c) * a.cross(b) + b.dot(b) * c.cross(a) +
              a.dot(a) * b.cross(c)) /
             det_2;
    radius = o.norm() + epsilon;
    center = O + o;
  }
  template <class Ary>
  Sphere(Ary const& pts) {
    *this = welzl_bounding_sphere(pts);
  }

  // Distance from p to boundary of the Sphere
  F signdis(Vec3 P) const { return (center - P).norm() - radius; }
  F signdis2(Vec3 P) const {  // NOT square of signdis!
    return (center - P).squaredNorm() - radius * radius;
  }

  bool intersect(Sphere other) const {
    F rtot = radius + other.radius;
    return (center - other.center).squaredNorm() <= rtot;
  }
  bool contact(Sphere other, F contact_dis) const {
    F rtot = radius + other.radius + contact_dis;
    return (center - other.center).squaredNorm() <= rtot * rtot;
  }
  bool contains(Vec3 pt) const {
    return (center - pt).squaredNorm() < radius * radius;
  }
};

template <class F>
constexpr F Sphere<F>::epsilon;

template <class F>
std::ostream& operator<<(std::ostream& out, Sphere<F> const& s) {
  out << s.center << " " << s.radius;
  return out;
}

template <class Ary>
auto welzl_bounding_sphere_impl(Ary const& points, size_t index,
                                std::vector<typename Ary::value_type>& sos,
                                size_t numsos) noexcept {
  using Pt = typename Ary::value_type;
  using F = typename Pt::Scalar;
  using Sph = Sphere<F>;
  // if no input points, the recursion has bottomed out. Now compute an
  // exact sphere based on points in set of support (zero through four points)
  if (index == 0) {
    switch (numsos) {
      case 0:
        return Sph();
      case 1:
        return Sph(sos[0]);
      case 2:
        return Sph(sos[0], sos[1]);
      case 3:
        return Sph(sos[0], sos[1], sos[2]);
      case 4:
        return Sph(sos[0], sos[1], sos[2], sos[3]);
    }
  }
  // Pick a point at "random" (here just the last point of the input set)
  --index;
  // Recursively compute the smallest bounding sphere of the remaining points
  Sph smallestSphere =
      welzl_bounding_sphere_impl(points, index, sos, numsos);  // (*)
  // If the selected point lies inside this sphere, it is indeed the smallest
  if (smallestSphere.contains(points[index])) return smallestSphere;
  // Otherwise, update set of support to additionally contain the new point
  assert(numsos < 4);
  sos[numsos] = points[index];
  // Recursively compute the smallest sphere of remaining points with new s.o.s.
  return welzl_bounding_sphere_impl(points, index, sos, numsos + 1);
}

template <class Ary>
auto welzl_bounding_sphere(Ary const& points) noexcept {
  using Pt = typename Ary::value_type;
  using F = typename Pt::Scalar;
  using Sph = Sphere<F>;
  std::vector<Pt> sos(4);
  return welzl_bounding_sphere_impl(points, points.size(), sos, 0);
}

template <class Ary>
auto central_bounding_sphere(Ary const& points) noexcept {
  using Pt = typename Ary::value_type;
  using F = typename Pt::Scalar;
  using Sph = Sphere<F>;
  Pt center;
  F radius = -1;
  if (points.size() > 0) {
    center = Pt(0, 0, 0);
    for (size_t i = 0; i < points.size(); i++) center += points[i];
    center /= (F)points.size();

    for (size_t i = 0; i < points.size(); i++) {
      F d2 = (points[i] - center).squaredNorm();
      if (d2 > radius) radius = d2;
    }
    radius = sqrt(radius) + Sph::epsilon;
  }
  return Sph(center, radius);
}
}
}
