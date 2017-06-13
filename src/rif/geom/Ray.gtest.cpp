#include <gtest/gtest.h>

#include <geom/Ray.hpp>

namespace rif {
namespace geom {
namespace test_numeric_Ray {

using V = V3<float>;
using X = X3<float>;

TEST(Ray, ray_basic_tests) {
  Ray<> r;
  // std::cout << r._m42 << std::endl;
  ASSERT_TRUE(r.orig().isApprox(V(0, 0, 0)));
  // std::cout << r.dirn() << std::endl;
  ASSERT_TRUE(r.dirn().isApprox(V(1, 0, 0)));
  ASSERT_FLOAT_EQ(r.dirn().norm(), 1.0);
  Ray<> s(V(5, 6, 7), V(2, 0, 0));
  ASSERT_TRUE(s.orig().isApprox(V(5, 6, 7)));
  ASSERT_TRUE(s.dirn().isApprox(V(1, 0, 0)));
  X trans10x = xform(V(10, 0, 0));
  auto t = trans10x * s;
  ASSERT_TRUE(t.dirn().isApprox(V(1, 0, 0)));
  ASSERT_TRUE(t.orig().isApprox(V(15, 6, 7)));
  t.orig() = V(4, 3, 2);
  t.dirn() = V(0, 0, 1);
  ASSERT_TRUE(t.orig().isApprox(V(4, 3, 2)));
  ASSERT_TRUE(t.dirn().isApprox(V(0, 0, 1)));
  X roty90(Eigen::AngleAxisf(0.5 * M_PI, V(0, 1, 0)));
  auto u = roty90 * t;
  ASSERT_TRUE(u.orig().isApprox(V(2, 3, -4)));
  ASSERT_TRUE(u.dirn().isApprox(V(1, 0, 0)));
  auto v = (trans10x * roty90) * u;
  ASSERT_TRUE(v.orig().isApprox(V(10 - 4, 3, -2)));
  ASSERT_TRUE(v.dirn().isApprox(V(0, 0, -1)));
}

TEST(Ray, rand_ray_gaussian) {
  for (int i = 0; i < 1000; ++i) {
    auto ray = rand_ray_gaussian<float>();
    ASSERT_FLOAT_EQ(ray.dirn().norm(), 1.0);
    ASSERT_LT(ray.orig().norm(), 200.0);
  }
}
}
}
}
