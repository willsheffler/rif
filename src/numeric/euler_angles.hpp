#ifndef INCLUDED_scheme_numeric_euler_angles_HH
#define INCLUDED_scheme_numeric_euler_angles_HH

#include <boost/math/constants/constants.hpp>
#include <cmath>

namespace scheme {
namespace numeric {

template <class M>
typename M::Scalar const &xx(M const &m) {
  return m(0, 0);
}
template <class M>
typename M::Scalar const &xy(M const &m) {
  return m(0, 1);
}
template <class M>
typename M::Scalar const &xz(M const &m) {
  return m(0, 2);
}
template <class M>
typename M::Scalar const &yx(M const &m) {
  return m(1, 0);
}
template <class M>
typename M::Scalar const &yy(M const &m) {
  return m(1, 1);
}
template <class M>
typename M::Scalar const &yz(M const &m) {
  return m(1, 2);
}
template <class M>
typename M::Scalar const &zx(M const &m) {
  return m(2, 0);
}
template <class M>
typename M::Scalar const &zy(M const &m) {
  return m(2, 1);
}
template <class M>
typename M::Scalar const &zz(M const &m) {
  return m(2, 2);
}
template <class M>
typename M::Scalar &xx(M &m) {
  return m(0, 0);
}
template <class M>
typename M::Scalar &xy(M &m) {
  return m(0, 1);
}
template <class M>
typename M::Scalar &xz(M &m) {
  return m(0, 2);
}
template <class M>
typename M::Scalar &yx(M &m) {
  return m(1, 0);
}
template <class M>
typename M::Scalar &yy(M &m) {
  return m(1, 1);
}
template <class M>
typename M::Scalar &yz(M &m) {
  return m(1, 2);
}
template <class M>
typename M::Scalar &zx(M &m) {
  return m(2, 0);
}
template <class M>
typename M::Scalar &zy(M &m) {
  return m(2, 1);
}
template <class M>
typename M::Scalar &zz(M &m) {
  return m(2, 2);
}

template <class M>
struct get_scalar {
  typedef typename M::Scalar type;
};

template <class T>
T sin_cos_range(T const &t) {
  return t;
}

/// @brief COPIED FROM ROSETTA
/// Return the three euler angles (in radians) that describe this
/// HomogeneousTransform as the series
/// of a Z axis rotation by the angle phi (returned in position 1 of the output
/// vector), followed by
/// an X axis rotation by the angle theta (returned in position 3 of the output
/// vector), followed by another
/// Z axis rotation by the angle psi (returned in position 2 of the output
/// vector).
/// This code is a modified version of Alex Z's code from r++.
///
/// @details
/// The range of phi is [ -pi, pi ];
/// The range of psi is [ -pi, pi ];
/// The range of theta is [ 0, pi ];
///
/// The function pretends that this HomogeneousTransform is the result of these
/// three transformations;
/// if it were, then the rotation matrix would be
///
/// FIGURE 1:
/// R = [
///       cos(psi)cos(phi)-cos(theta)sin(phi)sin(psi)
///       cos(psi)sin(phi)+cos(theta)cos(phi)sin(psi)      sin(psi)sin(theta)
///      -sin(psi)cos(phi)-cos(theta)sin(phi)cos(psi)
///      -sin(psi)sin(phi)+cos(theta)cos(phi)cos(psi)      cos(psi)sin(theta)
///                   sin(theta)sin(phi)
///                   -sin(theta)cos(phi)                        cos(theta)
/// ]
///
/// where each axis above is represented as a ROW VECTOR (to be distinguished
/// from the
/// HomogeneousTransform's representation of axes as COLUMN VECTORS).
///
/// The zz_ coordinate gives away theta.
/// Theta may be computed as acos( zz_ ), or, as Alex does it, asin( sqrt( 1 -
/// zz^2))
/// Since there is redundancy in theta, this function chooses a theta with a
/// positive
/// sin(theta): i.e. quadrants I and II.  Assuming we have a positive sin theta
/// pushes phi and psi into conforming angles.
///
/// NOTE on theta: asin returns a value in the range [ -pi/2, pi/2 ], and we
/// have artificially
/// created a positive sin(theta), so we will get a asin( pos_sin_theta ), we
/// have a value
/// in the range [ 0, pi/2 ].  To convert this into the actual angle theta, we
/// examine the zz sign.
/// If zz is negative, we chose the quadrant II theta.
/// That is, asin( pos_sin_theta) returned an angle, call it theta'.  Now, if
/// cos( theta ) is negative,
/// then we want to choose the positive x-axis rotation that's equivalent to
/// -theta'.  To do so,
/// we reflect q through the y axis (See figure 2 below) to get p and then
/// measure theta as pi - theta'.
///
/// FIGURE 2:
///
///  II        |         I
///            |
///    p.      |      .q (cos(-theta'), abs(sin(theta')))
///       .    |    .
/// theta'( .  |  .  )  theta' = asin( abs(sin(theta))
/// -----------------------
///            |
///            |
///            |
///  III       |        IV
///            |
///  The angle between the positive x axis and p is pi - theta'.
///
///
///
/// Since zx and zy contain only phi terms and a constant sin( theta ) term,
/// phi is given by atan2( sin_phi, cos_phi ) = atan2( c*sin_phi, c*cos_phi ) =
/// atan2( zx, -zy )
/// for c positive and non-zero.  If sin_theta is zero, or very close to zero,
/// we're at gimbal lock.
///
/// Moreover, since xz and yz contain only psi terms, psi may also be deduced
/// using atan2.
///
/// There are 2 degenerate cases (gimbal lock)
/// 1. theta close to 0  (North Pole singularity), or
/// 2. theta close to pi (South Pole singularity)
/// For these, we take: phi=acos(xx), theta = 0 (resp. Pi/2), psi = 0
template <class M, class E>
void euler_angles(M const &m, E &euler) {
  typedef typename get_scalar<M>::type T;
  static T const pi = boost::math::constants::pi<T>();
  static T const pi_2 = boost::math::constants::pi<T>() * (T)2;
  static T const epsilon = (T)10 * std::numeric_limits<T>::epsilon();
  static T const FLOAT_PRECISION = std::sqrt(epsilon);
  if (zz(m) >= (T)1 - FLOAT_PRECISION) {
    euler[0] = std::atan2(sin_cos_range(yx(m)), sin_cos_range(xx(m)));
    euler[1] = 0.0;
    euler[2] = 0.0;
  } else if (zz(m) <= (T)-1 + FLOAT_PRECISION) {
    euler[0] = std::atan2(sin_cos_range(yx(m)), sin_cos_range(xx(m)));
    euler[1] = 0.0;
    euler[2] = M_PI;
  } else {
    T pos_sin_theta =
        std::sqrt((T)1 - zz(m) * zz(m));  // sin2theta = 1 - cos2theta.
    euler[2] = std::asin(pos_sin_theta);
    if (zz(m) < 0) {
      euler[2] = pi - euler[2];
    }
    euler[0] = std::atan2(xz(m), -yz(m));
    euler[1] = std::atan2(zx(m), zy(m));
  }
  euler[0] += euler[0] < 0.0 ? pi_2 : 0.0;
  euler[1] += euler[1] < 0.0 ? pi_2 : 0.0;

  euler[0] = std::min(std::max(0.0, euler[0]), pi_2 - epsilon);
  euler[1] = std::min(std::max(0.0, euler[1]), pi_2 - epsilon);
  euler[2] = std::min(std::max(0.0, euler[2]), pi - epsilon);

  assert(0 <= euler[0]);
  assert(euler[0] < pi_2);
  assert(0 <= euler[1]);
  assert(euler[1] < pi_2);
  assert(0 <= euler[2]);
  assert(euler[2] <= pi);
}

///@brief euler to matrix
template <class M, class E>
void from_euler_angles(E const &euler, M &m) {
  typedef typename get_scalar<M>::type T;
  T const ce1(std::cos(euler[0])), se1(std::sin(euler[0]));
  T const ce2(std::cos(euler[1])), se2(std::sin(euler[1]));
  T const ce3(std::cos(euler[2])), se3(std::sin(euler[2]));
  xx(m) = ce2 * ce1 - ce3 * se1 * se2;
  yx(m) = ce2 * se1 + ce3 * ce1 * se2;
  zx(m) = se2 * se3;
  xy(m) = -se2 * ce1 - ce3 * se1 * ce2;
  yy(m) = -se2 * se1 + ce3 * ce1 * ce2;
  zy(m) = ce2 * se3;
  xz(m) = se3 * se1;
  yz(m) = -se3 * ce1;
  zz(m) = ce3;
}

template <class M, class E>
void euler_angles_deg(M const &m, E &euler) {
  typedef typename get_scalar<M>::type T;
  static T const rad_to_deg = 180.0 / boost::math::constants::pi<T>();
  euler_angles(m, euler);
  euler[0] = euler[0] * rad_to_deg;
  euler[1] = euler[1] * rad_to_deg;
  euler[2] = euler[2] * rad_to_deg;
}
template <class M, class E>
void from_euler_angles_deg(E const &erad, M &m) {
  typedef typename get_scalar<M>::type T;
  static T const deg_to_rad = boost::math::constants::pi<T>() / 180.0;
  E euler;
  euler[0];
  euler[0] = erad[0] * deg_to_rad;
  euler[1] = erad[1] * deg_to_rad;
  euler[2] = erad[2] * deg_to_rad;
  euler_angles(m, euler);
}
}
}

#endif
