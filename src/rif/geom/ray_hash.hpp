#pragma once

#include <eigen_types.hpp>
#include <geom/Ray.hpp>
#include <geom/cube_to_sphere.hpp>
#include <iostream>
#include <numeric/bcc_lattice.hpp>
#include <util/SimpleArray.hpp>

namespace rif {
namespace geom {

using namespace rif::numeric;

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
 * @detail     canonical position places r at the orig along x and s in the
 * x,y plane
 */
template <class F>
Ray<F> align_ray_pair(Ray<F> a, Ray<F> b) {
  // todo: should this be made symmetrical???
  M3<F> rotation;
  V3<F> basis1 = a.dirn;
  V3<F> basis3 = basis1.cross(b.orig - a.orig).normalized();
  V3<F> basis2 = basis3.cross(basis1).normalized();
  rotation.row(0) = basis1;
  rotation.row(1) = basis2;
  rotation.row(2) = basis3;
  V3<F> a2_origin = rotation * a.orig;
  Ray<F> b2 = rotation * b;
  b2.orig -= a2_origin;
  // std::cout << rotation.determinant() << std::endl;
  // Ray<F> a2 = rotation * a;
  // std::cout << "a  " << a << std::endl;
  // std::cout << "b  " << b << std::endl;
  // std::cout << "a2 " << a2 << std::endl;
  // std::cout << "b2 " << b2 << std::endl;
  return b2;
}

///@brief get num cells for quadsphere
template <class F>
int qs_nc(F resl, F lever) {
  // TODO: audit this fudge factor
  return std::max(F(1), lever / resl * F(1.5));
}
template <class F>
F qs_bound(F resl, F lever) {
  return 1.0 - 0.25 / qs_nc(resl, lever);
}

template <class RHash>
auto brute_maxmin_nbr(RHash& rh, bool face0_only = false) {
  using R = typename RHash::R;
  using F = typename RHash::F;
  int nface = face0_only ? 1 : 6;
  F maxmindot = -9e9, maxmindis = -9e9;
  for (int i = 0; i < rh.bcc_[0].size() * nface - 1; ++i) {
    R a = rh.get_center(i);
    F mindot = 9e9;
    F mindis = 9e9;
    for (int j = i + 1; j < rh.bcc_[0].size() * nface; ++j) {
      R b = rh.get_center(j);
      mindot = std::min(mindot, a.dirn.dot(b.dirn));
      F dis = (a.orig - b.orig).norm();
      if (dis > 0.0001) mindis = std::min(mindis, dis);
    }
    maxmindot = std::max(maxmindot, mindot);
    maxmindis = std::max(maxmindis, mindis);
  }
  return A2<F>(maxmindis, acos(maxmindot));
}

/**
 * @brief      Bin relation between of rays (4D)
 *
 * @tparam     _R    { description }
 * @tparam     _K    { description }
 * @detail align_ray_pair to put ray2 in 'canonical' position
 */
template <class _R = Ray<float>, class _K = uint32_t>
struct RayToRay4dHash {
  using R = _R;
  using F = typename R::Scalar;
  using K = _K;
  using Grid = BCC<4, F, K>;
  using Fs = typename Grid::Floats;
  using V = typename R::V;
  util::SimpleArray<6, Grid> bcc_;
  F resl_, lever_, bound_;
  RayToRay4dHash(F resl = 0.25, F lever = 1.5, F bound = 128)
      : bcc_(Grid(
            I4(2 * bound / resl, 2 * bound / resl, qs_nc(resl, lever),
               qs_nc(resl, lever)),
            Fs(-bound, -bound, -qs_bound(resl, lever), -qs_bound(resl, lever)),
            Fs(+bound, +bound, +qs_bound(resl, lever),
               +qs_bound(resl, lever)))),
        resl_(resl),
        lever_(lever),
        bound_(bound) {}
  K get_key(R a, R b) const noexcept { return get_key(align_ray_pair(a, b)); }
  K get_key(R r) const noexcept {
    assert(fabs(r.orig[2]) < 0.001);
    // std::cout << "get_key: " << r << std::endl;
    F qsA, qsB;
    K face = get_quadsphere_coords(r.dirn, qsA, qsB);
    Fs bcrd(r.orig[0], r.orig[1], qsA, qsB);
    // std::cout << "get_key: " << bcrd << std::endl;
    K bcc_key = bcc_[face][bcrd];
    // std::cout << "get_key: " << face << " " << bcc_key << std::endl;
    return face * bcc_[face].size() + bcc_key;
  }
  R get_center(K k) const noexcept {
    K face = k / bcc_[face].size();
    K bcc_key = k % bcc_[face].size();
    // std::cout << "get_cen: " << face << " " << bcc_key << std::endl;
    Fs bcrd = bcc_[face][bcc_key];
    // std::cout << "get_cen: " << bcrd << std::endl;
    auto dir = get_point_from_quadsphere_coords<V>(face, bcrd[2], bcrd[3]);
    R r(V(bcrd[0], bcrd[1], 0.0), dir);
    return r;
  }
  uint size() const noexcept { return 6 * bcc_[0].size(); }
  uint size_cart() const noexcept { return bcc_[0].nside_[0]; }
  uint size_qsph() const noexcept { return bcc_[0].nside_[2]; }
};

/**
 * @brief      Bin Rays 5D
 *
 * @tparam     _R    { description }
 * @tparam     _K    { description }
 */
template <class _R = Ray<float>, class _K = uint32_t>
struct Ray5dHash {
  // todo: reduce code duplication with RayToRay4dHash
  using R = _R;
  using F = typename R::Scalar;
  using K = _K;
  using Grid = BCC<5, F, K>;
  using Fs = typename Grid::Floats;
  using V = typename R::V;
  util::SimpleArray<6, Grid> bcc_;
  F resl_, lever_, bound_;
  Ray5dHash(F resl = 0.25, F lever = 1.5, F bound = 32)
      : bcc_(Grid(makeI5(2 * bound / resl, 2 * bound / resl, 2 * bound / resl,
                         qs_nc(resl, lever), qs_nc(resl, lever)),
                  Fs(-bound, -bound, -bound, -qs_bound(resl, lever),
                     -qs_bound(resl, lever)),
                  Fs(+bound, +bound, +bound, +qs_bound(resl, lever),
                     +qs_bound(resl, lever)))),
        resl_(resl),
        lever_(lever),
        bound_(bound) {}
  K get_key(R r) const noexcept {
    // std::cout << "get_key: " << r << std::endl;
    F qsA, qsB;
    K face = get_quadsphere_coords(r.dirn, qsA, qsB);
    Fs bcrd(r.orig[0], r.orig[1], r.orig[2], qsA, qsB);
    // std::cout << "get_key: " << bcrd << std::endl;
    K bcc_key = bcc_[face][bcrd];
    // std::cout << "get_key: " << face << " " << bcc_key << std::endl;
    return face * bcc_[face].size() + bcc_key;
  }
  R get_center(K k) const noexcept {
    K face = k / bcc_[0].size();
    K bcc_key = k % bcc_[0].size();
    // std::cout << "get_cen: " << face << " " << bcc_key << std::endl;
    auto bcrd = bcc_[face][bcc_key];
    // std::cout << "get_cen: " << bcrd << std::endl;
    auto dir = get_point_from_quadsphere_coords<V>(face, bcrd[3], bcrd[4]);
    V v;
    R r(V(bcrd[0], bcrd[1], bcrd[2]), dir);
    return r;
  }
  uint size() const noexcept { return 6 * bcc_[0].size(); }
  uint size_cart() const noexcept { return bcc_[0].nside_[0]; }
  uint size_qsph() const noexcept { return bcc_[0].nside_[3]; }
};

/**
 * @brief      Bin Pairs of Rays in space 10D
 *
 * @tparam     _R    { description }
 * @tparam     _K    { description }
 */
template <class _R = Ray<float>, class _K = uint64_t>
struct RayRay10dHash {
  // todo: reduce code duplication with RayToRay4dHash
  using R = _R;
  using F = typename R::Scalar;
  using K = _K;
  using Grid = BCC<10, F, K>;
  using Fs = typename Grid::Floats;
  using V = typename R::V;
  util::SimpleArray<36, Grid> bcc_;
  F resl_, lever_, bound_;
  RayRay10dHash(F resl = 0.25, F lever = 1.5, F bound = 32)
      : bcc_(Grid(makeI10(2 * bound / resl, 2 * bound / resl, 2 * bound / resl,
                          qs_nc(resl, lever), qs_nc(resl, lever),
                          2 * bound / resl, 2 * bound / resl, 2 * bound / resl,
                          qs_nc(resl, lever), qs_nc(resl, lever)),
                  Fs(-bound, -bound, -bound, -qs_bound(resl, lever),
                     -qs_bound(resl, lever), -bound, -bound, -bound,
                     -qs_bound(resl, lever), -qs_bound(resl, lever)),
                  Fs(+bound, +bound, +bound, +qs_bound(resl, lever),
                     +qs_bound(resl, lever), +bound, +bound, +bound,
                     +qs_bound(resl, lever), +qs_bound(resl, lever)))),
        resl_(resl),
        lever_(lever),
        bound_(bound) {}
  K get_key(R r, R s) const noexcept {
    // std::cout << "get_key: " << r << std::endl;
    F qsx1, qsy1;
    F qsx2, qsy2;
    K face1 = get_quadsphere_coords(r.dirn, qsx1, qsy1);
    K face2 = get_quadsphere_coords(s.dirn, qsx2, qsy2);
    K face12 = 6 * face1 + face2;
    Fs bcrd(r.orig[0], r.orig[1], r.orig[2], qsx1, qsy1, s.orig[0], s.orig[1],
            s.orig[2], qsx2, qsy2);
    // std::cout << "get_key: " << bcrd << std::endl;
    K bcc_key = bcc_[face12][bcrd];
    // std::cout << "get_key: " << face << " " << bcc_key << std::endl;
    return face12 * bcc_[0].size() + bcc_key;
  }
  K get_key_from_pair(std::pair<R, R> p) const noexcept {
    return get_key(p.first, p.second);
  }
  std::pair<R, R> get_center(K k) const {
    K face12 = k / bcc_[0].size();
    K face1 = face12 / 6;  // !!!!!!!!!!!!! 1/2 swap
    K face2 = face12 % 6;
    K bcc_key = k % bcc_[0].size();
    // std::cout << "get_cen: " << face << " " << bcc_key << std::endl;
    auto bcrd = bcc_[face12][bcc_key];
    // std::cout << "get_cen: " << bcrd << std::endl;
    auto dir1 = get_point_from_quadsphere_coords<V>(face1, bcrd[3], bcrd[4]);
    auto dir2 = get_point_from_quadsphere_coords<V>(face2, bcrd[8], bcrd[9]);
    V v;
    R r(V(bcrd[0], bcrd[1], bcrd[2]), dir1);
    R s(V(bcrd[5], bcrd[6], bcrd[7]), dir2);
    return std::make_pair(r, s);
  }
  uint size() const noexcept { return 36 * bcc_[0].size(); }
  uint size_cart() const noexcept { return bcc_[0].nside_[0]; }
  uint size_qsph() const noexcept { return bcc_[0].nside_[3]; }
};
}
}
