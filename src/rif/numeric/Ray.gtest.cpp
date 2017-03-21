#include <gtest/gtest.h>

#include <numeric/Ray.hpp>

namespace rif {
namespace numeric {
namespace test_numeric_Ray {

TEST(Ray, Ray_basic_tests) {
  Ray r;
  ASSERT_FLOAT_EQ(r.direction.norm(), 1.0);
  Ray s(V(0, 0, 0), V(2, 0, 0));
  ASSERT_FLOAT_EQ(s.direction.norm(), 1.0);
  X trans10 = xform(V(10, 0, 0));
  Ray t = trans10 * s;
  ASSERT_TRUE(t.direction.isApprox(V(1, 0, 0)));
  ASSERT_TRUE(t.origin.isApprox(V(10, 0, 0)));
}
}
}
}
