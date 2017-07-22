#include <gtest/gtest.h>

#include "rif/geom/ray_hash.hpp"
#include "rif/numeric/global_rng.hpp"

namespace rif {
namespace geom {
namespace ray_hash_test {

using std::cout;
using std::endl;
using rif::numeric::global_rng;

TEST(ray_hash, align_ray_pair) {
  // std::mt19937 rng((unsigned int)time(0) + 267943);
  for (int i = 0; i < 100; ++i) {
    auto a1 = rand_ray_gaussian(/*rng*/);
    auto a2 = rand_ray_gaussian(/*rng*/);
    auto b1 = Ray<>(V3f(0, 0, 0), V3f(1, 0, 0));
    auto b2 = align_ray_pair(a1, a2);
    ASSERT_LT(fabs(b2.orig()[2]), 0.001f);
    ASSERT_NEAR(a1.dirn().dot(a2.dirn()), b1.dirn().dot(b2.dirn()), 0.0001);
    ASSERT_NEAR((a1.orig() - a2.orig()).norm(), (b1.orig() - b2.orig()).norm(),
                0.0001);
    ASSERT_NEAR((a1.orig() + a1.dirn() - a2.orig() - a2.dirn()).norm(),
                (b1.orig() + b1.dirn() - b2.orig() - b2.dirn()).norm(), 0.0001);
  }
}

TEST(ray_hash, RayToRay4dHash_invertibility) {
  auto rh = RayToRay4dHash<>(0.25, 1.0, 1.0);
  double nerr = 0;
  for (int key = 0; key < rh.size(); ++key) {
    auto cen = rh.get_center(key);
    auto cen_key = rh.get_key_aligned(cen);
    if (key != cen_key) {
      if (rh.quadsphere_inbounds(key)) ASSERT_EQ(key, cen_key);
    }
  }
}

TEST(ray_hash, Ray5dHash_invertibility) {
  Ray5dHash<> rh(0.33, 1.0, 1.0);
  // cout << "qsph ncell " << rh.size_qsph() << endl;
  double nerr = 0;
  for (int key = 0; key < rh.size(); ++key) {
    auto cen = rh.get_center(key);
    auto cen_key = rh.get_key(cen);
    if (key != cen_key) {
      if (rh.quadsphere_inbounds(key)) ASSERT_EQ(key, cen_key);
    }
  }
}

TEST(ray_hash, RayRay10dHash_invertibility) {
  RayRay10dHash<> rh(0.33, 1.0, 1.0);
  double nerr = 0;
  // std::cout << rh.size() << std::endl;
  uint64_t step = sqrt(rh.size()) / 16;
  // cout << "nsamp " << rh.size() / step << endl;
  for (uint64_t key = 0; key < rh.size(); key += step) {
    auto cen = rh.get_center(key);
    auto cen_key = rh.get_key_from_pair(cen);
    nerr += key != cen_key;
    if (key != cen_key) {
      if (rh.quadsphere_inbounds(key)) ASSERT_EQ(key, cen_key);
    }
  }
}

TEST(ray_hash, ray_ray_hash_cart_ori_spacing) {
  RayToRay4dHash<> rh(0.25, 1.0, 0.3);
  A2f nbr_sp = brute_maxmin_nbr(rh, 1);
  // std::cout << "rh.size " << rh.size() << ", cart " << rh.size_cart()
  // << ", qsph " << rh.size_qsph()
  // << ", closest: " << nbr_sp.transpose() << std::endl;
  ASSERT_GE(nbr_sp[0], 1.5 * nbr_sp[1]);
}

TEST(ray_hash, ray_hash_cart_ori_spacing) {
  Ray5dHash<> rh(0.4, 1.0, 0.5);
  A2f nbr_sp = brute_maxmin_nbr(rh, 1);
  std::cout << "             rh.size " << rh.size() << ", cart "
            << rh.size_cart() << ", qsph " << rh.size_qsph()
            << ", closest: " << nbr_sp.transpose() << std::endl;
  ASSERT_GE(nbr_sp[0], 1.4 * nbr_sp[1]);
}
}
}
}
