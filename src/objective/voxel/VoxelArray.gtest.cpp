#include <gtest/gtest.h>

#include "io/cache.hpp"
#include "objective/voxel/VoxelArray.hpp"

#include <boost/foreach.hpp>
#include <random>

namespace scheme {
namespace objective {
namespace voxel {
namespace test {

using std::cout;
using std::endl;

TEST(VoxelArray, simple1d) {
  typedef util::SimpleArray<1, float> F1;
  // VoxelArray<0,float> a; // this should fail to compile

  VoxelArray<1, float> a1(0, 3, 0.5);
  ASSERT_EQ(a1.size(), 7);
  a1[F1(0)] = 1;
  ASSERT_EQ(a1[F1(0)], 1);
#ifdef DEBUG
  // ASSERT_DEATH( a1[F1(-0.1)], ".*" );
  a1[F1(0)];
  a1[F1(2.999)];
// ASSERT_DEATH( a1[F1(3.0)], ".*" );
#endif
}
TEST(VoxelArray, simple3d) {
  typedef util::SimpleArray<3, float> F3;
  VoxelArray<3, float> a3(F3(-1, -2, -3), F3(1, 2, 3), 0.5);
  ASSERT_EQ(a3.shape()[0], 5);
  ASSERT_EQ(a3.shape()[1], 9);
  ASSERT_EQ(a3.shape()[2], 13);
  a3[F3(0, 0, 0)] = 1;
  a3[F3(-1, 1, 2.3)] = 2;
  ASSERT_EQ(a3[F3(0, 0, 0)], 1);
#ifdef DEBUG
// ASSERT_DEATH( a3[F3(1,2,3)], ".*" );
// ASSERT_DEATH( a3[F3(-1,-2,-3.1)], ".*" );
#endif
  ASSERT_EQ(a3[F3(0, 0, 0)], 1);
  ASSERT_EQ(a3[F3(-1, 1, 2.3)], 2);
}
TEST(VoxelArray, bounds3d) {
  typedef util::SimpleArray<3, float> F3;
  VoxelArray<3, float> a3(F3(-1, -2, -3), F3(1, 2, 3), 0.49);
  ASSERT_EQ(a3.shape()[0], 5);
  ASSERT_EQ(a3.shape()[1], 9);
  ASSERT_EQ(a3.shape()[2], 13);
  a3[F3(0, 0, 0)] = 1;
  a3[F3(-1, 1, 2.3)] = 2;
  ASSERT_EQ(a3[F3(0, 0, 0)], 1);
#ifdef DEBUG
// ASSERT_DEATH( a3[F3(1,2,3)], ".*" );
// ASSERT_DEATH( a3[F3(-1,-2,-3.1)], ".*" );
#endif
  ASSERT_EQ(a3[F3(0, 0, 0)], 1);
  ASSERT_EQ(a3[F3(-1, 1, 2.3)], 2);
}

TEST(VoxelArray, io) {
  std::mt19937 rng((unsigned int)time(0));
  std::uniform_real_distribution<> uniform;

  VoxelArray<3, double> a(-6, 7, 1.6345);
  for (size_t i = 0; i < a.num_elements(); ++i) a.data()[i] = uniform(rng);

  VoxelArray<3, double> b(-6, 7, 1.6354);
  ASSERT_NE(a, b);
  b = a;
  ASSERT_EQ(a, b);
  for (size_t i = 0; i < b.num_elements(); ++i) b.data()[i] = 0;
  ASSERT_NE(a, b);

#ifdef CEREAL
  ASSERT_EQ(a, io::test_serialization(a));
#endif

  std::ostringstream oss;
  a.save(oss);
  std::istringstream iss(oss.str());
  VoxelArray<3, double> save_read_a;
  save_read_a.load(iss);

  ASSERT_TRUE(a == save_read_a);
}
}
}
}
}
