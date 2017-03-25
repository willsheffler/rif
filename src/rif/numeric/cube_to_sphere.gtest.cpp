#include <gtest/gtest.h>

#include <random>

#include "numeric/cube_to_sphere.hpp"
#include "util/SimpleArray.hpp"

namespace rif {
namespace numeric {
namespace cube_to_sphere_test {

using V3 = util::SimpleArrayLegacy<3, float>;
std::mt19937 rng((unsigned int)time(0) + 2746387);
std::normal_distribution<> rnorm;
using std::cout;
using std::endl;

TEST(cube_to_sphere, test_permute_and_inv) {
  std::uniform_int_distribution<> randint(0, 5);
  for (int i = 0; i < 1000; ++i) {
    int facenum = randint(rng);
    float x, x0, y, y0, z, z0;
    x0 = x = rnorm(rng);
    y0 = y = rnorm(rng);
    z0 = z = rnorm(rng);
    permute_cube_face(facenum, x, y, z);
    permute_cube_face_inverse(facenum, x, y, z);
    ASSERT_FLOAT_EQ(x0, x);
    ASSERT_FLOAT_EQ(y0, y);
    ASSERT_FLOAT_EQ(z0, z);

    V3 v0(x, y, z);
    v0.normalize();
    V3 v = v0;
    permute_cube_face(facenum, v);
    permute_cube_face_inverse(facenum, v);
    ASSERT_TRUE(v0.isApprox(v));
  }
}
// order supposed to be Z,X,-Y,-X,Y,-Z
TEST(cube_to_sphere, get_cube_facenum) {
  ASSERT_EQ(0, get_cube_facenum(V3(+0, +0, +1)));
  ASSERT_EQ(1, get_cube_facenum(V3(+1, +0, +0)));
  ASSERT_EQ(2, get_cube_facenum(V3(+0, -1, +0)));
  ASSERT_EQ(3, get_cube_facenum(V3(-1, +0, +0)));
  ASSERT_EQ(4, get_cube_facenum(V3(+0, +1, +0)));
  ASSERT_EQ(5, get_cube_facenum(V3(+0, +0, -1)));
}

TEST(cube_to_sphere, cube_to_sphere_and_back) {
  for (int i = 0; i < 1000; ++i) {
    V3 v(rnorm(rng), rnorm(rng), rnorm(rng));
    v.normalize();
    auto v0 = v;
    auto face = get_cube_facenum(v);
    ASSERT_LT(face, 6);
    ASSERT_GT(-1, face);
    // cout << face << " " << v << endl;
    permute_cube_face_inverse(face, v);
    ASSERT_EQ(get_cube_facenum(v), 0);
    auto v1 = v;
    sphere_to_cube_facenum0(v);
    auto v2 = v;
    // cout << 0 << " " << v << endl;
    ASSERT_FLOAT_EQ(v[2], 1.0);
    permute_cube_face(face, v);

    cube_to_sphere(v);
    if (!v0.isApprox(v))
      cout << "v0 " << v0 << "\nv1 " << v1 << "\nv2 " << v2 << "\nv  " << v
           << endl;
    ASSERT_TRUE(v0.isApprox(v));
  }
}
}
}
}