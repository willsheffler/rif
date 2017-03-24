#include <gtest/gtest.h>

#include "numeric/ray_ray_hash.hpp"

namespace rif {
namespace numeric {
namespace ray_ray_hash_test {

TEST(ray_ray_hash, align_ray_pair) {
  std::mt19937 rng((unsigned int)time(0) + 267943);
  for (int i = 0; i < 100; ++i) {
    auto a1 = rand_ray_gaussian(rng);
    auto a2 = rand_ray_gaussian(rng);
    auto b1 = Ray<>(V(0, 0, 0), V(1, 0, 0));
    auto b2 = align_ray_pair(a1, a2);
    ASSERT_LT(fabs(b2.origin[2]), 0.001f);
    ASSERT_NEAR(a1.direction.dot(a2.direction), b1.direction.dot(b2.direction),
                0.0001);
    ASSERT_NEAR((a1.origin - a2.origin).norm(), (b1.origin - b2.origin).norm(),
                0.0001);
  }
}

TEST(ray_ray_hash, basic_test) {
  std::cout << "WRITE THIS SHIT!!!" << std::endl;
  ASSERT_TRUE(1);
}
}
}
}
