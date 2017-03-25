#include <gtest/gtest.h>

#include <numeric/Ray.hpp>

namespace rif {
namespace numeric {
namespace test_numeric_Ray {

using V = V3<float>;
using X = X3<float>;

TEST(Ray, ray_basic_tests) {
  Ray<> r;
  ASSERT_FLOAT_EQ(r.direction.norm(), 1.0);
  Ray<> s(V(0, 0, 0), V(2, 0, 0));
  ASSERT_FLOAT_EQ(s.direction.norm(), 1.0);
  X trans10 = xform(V(10, 0, 0));
  Ray<> t = trans10 * s;
  ASSERT_TRUE(t.direction.isApprox(V(1, 0, 0)));
  ASSERT_TRUE(t.origin.isApprox(V(10, 0, 0)));
}

TEST(Ray, rand_ray_gaussian) {
  std::mt19937 rng((unsigned int)time(0) + 267943);
  for (int i = 0; i < 10; ++i) {
    auto ray = rand_ray_gaussian(rng);
    ASSERT_FLOAT_EQ(ray.direction.norm(), 1.0);
    ASSERT_LT(ray.origin.norm(), 50.0);
  }
}
}
}
}
