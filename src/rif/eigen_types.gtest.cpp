#include "eigen_types.hpp"
#include <gtest/gtest.h>

namespace rif {
namespace eigen_types_test {

TEST(eigen_types, utility_funcs) {
  // I1 i1(3);  // this doesn't work
  I1 i1 = makeI1(3);
  for (int i = 0; i < 1; ++i) ASSERT_EQ(i1[i], 3 + i);
  I2 i2(3, 4);
  for (int i = 0; i < 2; ++i) ASSERT_EQ(i2[i], 3 + i);
  I3 i3(3, 4, 5);
  for (int i = 0; i < 3; ++i) ASSERT_EQ(i3[i], 3 + i);
  I4 i4(3, 4, 5, 6);
  for (int i = 0; i < 4; ++i) ASSERT_EQ(i4[i], 3 + i);
  I5 i5 = makeI5(3, 4, 5, 6, 7);
  for (int i = 0; i < 5; ++i) ASSERT_EQ(i5[i], 3 + i);
  I6 i6 = makeI6(3, 4, 5, 6, 7, 8);
  for (int i = 0; i < 6; ++i) ASSERT_EQ(i6[i], 3 + i);
  I7 i7 = makeI7(3, 4, 5, 6, 7, 8, 9);
  for (int i = 0; i < 7; ++i) ASSERT_EQ(i7[i], 3 + i);
  I8 i8 = makeI8(3, 4, 5, 6, 7, 8, 9, 10);
  for (int i = 0; i < 8; ++i) ASSERT_EQ(i8[i], 3 + i);
  I9 i9 = makeI9(3, 4, 5, 6, 7, 8, 9, 10, 11);
  for (int i = 0; i < 9; ++i) ASSERT_EQ(i9[i], 3 + i);
  I10 i10 = makeI10(3, 4, 5, 6, 7, 8, 9, 10, 11, 12);
  for (int i = 0; i < 10; ++i) ASSERT_EQ(i10[i], 3 + i);
}
}
}