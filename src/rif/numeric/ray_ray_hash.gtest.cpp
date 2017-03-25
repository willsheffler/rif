#include <gtest/gtest.h>

#include "numeric/ray_ray_hash.hpp"

namespace rif {
namespace numeric {
namespace ray_ray_hash_test {

using std::cout;
using std::endl;

TEST(ray_ray_hash, align_ray_pair) {
  std::mt19937 rng((unsigned int)time(0) + 267943);
  for (int i = 0; i < 100; ++i) {
    auto a1 = rand_ray_gaussian(rng);
    auto a2 = rand_ray_gaussian(rng);
    auto b1 = Ray<>(V3f(0, 0, 0), V3f(1, 0, 0));
    auto b2 = align_ray_pair(a1, a2);
    ASSERT_LT(fabs(b2.origin[2]), 0.001f);
    ASSERT_NEAR(a1.direction.dot(a2.direction), b1.direction.dot(b2.direction),
                0.0001);
    ASSERT_NEAR((a1.origin - a2.origin).norm(), (b1.origin - b2.origin).norm(),
                0.0001);
    ASSERT_NEAR((a1.origin + a1.direction - a2.origin - a2.direction).norm(),
                (b1.origin + b1.direction - b2.origin - b2.direction).norm(),
                0.0001);
  }
}

TEST(ray_ray_hash, test_bins_of_centers) {
  RayBins<> rh(0.25, 1.0, 1.0);
  double nerr = 0;
  for (int key = 0; key < rh.size(); ++key) {
    auto cen = rh.get_center(key);
    auto cen_key = rh.get_key(cen);
    nerr += key != cen_key;
    ASSERT_EQ(key, cen_key);
  }
  ASSERT_LT(nerr / rh.size(), 0.1);
}

TEST(ray_ray_hash, test_cart_ori_spacing) {
  RayBins<> rh(0.25, 1.0, 0.3);
  A2<float> nbr_sp = rh.brute_maxmin_nbr(1);
  // std::cout << "rh.size " << rh.size() << ", cart " << rh.size_cart()
  // << ", qsph " << rh.size_qsph()
  // << ", closest: " << nbr_sp.transpose() << std::endl;
  ASSERT_GE(nbr_sp[0], 1.5 * nbr_sp[1]);
}
}
}
}
