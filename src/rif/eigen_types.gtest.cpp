#include <gtest/gtest.h>
#include "eigen_types.hpp"

namespace rif {
namespace eigen_types_test {

using std::cout;
using std::endl;

TEST(eigen_types, ctor_funcs) {
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

TEST(eigen_types, raw_asif_rowmajor) {
  using util::raw_asif_rowmajor;
  using namespace Eigen;
  using X4r = Transform<float, 3, Affine, RowMajor>;
  using X4c = Transform<float, 3, Affine, ColMajor>;
  using X3r = Transform<float, 3, AffineCompact, RowMajor>;
  using X3c = Transform<float, 3, AffineCompact, ColMajor>;
  auto x4r = X4r::Identity();
  auto x4c = X4c::Identity();
  auto x3r = X3c::Identity();
  auto x3c = X3c::Identity();

  M3f rot;
  rot << 10, 11, 12, 13, 14, 15, 16, 17, 18;

  x4r.linear() = rot;
  x4c.linear() = rot;
  x3r.linear() = rot;
  x3c.linear() = rot;

  x4r.translation() = V3f(20, 21, 22);
  x4c.translation() = V3f(20, 21, 22);
  x3r.translation() = V3f(20, 21, 22);
  x3c.translation() = V3f(20, 21, 22);

  cout << "=== X3r === (is ColMajor (in spite of asking for RowMajor)!!!!)"
       << endl;
  for (int i = 0; i < sizeof(X3r) / sizeof(float); ++i)
    cout << ' ' << x3r.data()[i];
  cout << endl;
  for (int i = 0; i < 4; ++i) {
    for (int j = 0; j < 3; ++j) {
      auto d = x3r.data()[i * 3 + j];
      cout << (d > 9 ? "" : " ") << d << ' ';
    }
    cout << endl;
  }
  cout << endl;

  cout << "=== X3c === " << endl;
  ;
  for (int i = 0; i < sizeof(X3c) / sizeof(float); ++i)
    cout << ' ' << x3c.data()[i];
  cout << endl;
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 4; ++j) {
      auto d = x3c.data()[i * 4 + j];
      cout << (d > 9 ? "" : " ") << d << ' ';
    }
    cout << endl;
  }
  cout << endl;

  cout << "=== X4r ===" << endl;
  for (int i = 0; i < sizeof(X4r) / sizeof(float); ++i)
    cout << ' ' << x4r.data()[i];
  cout << endl;
  for (int i = 0; i < 4; ++i) {
    for (int j = 0; j < 4; ++j) {
      auto d = x4r.data()[i * 4 + j];
      cout << (d > 9 ? "" : " ") << d << ' ';
    }
    cout << endl;
  }
  cout << endl;

  cout << "=== X4c ===" << endl;
  for (int i = 0; i < sizeof(X4c) / sizeof(float); ++i)
    cout << ' ' << x4c.data()[i];
  cout << endl;
  for (int i = 0; i < 4; ++i) {
    for (int j = 0; j < 4; ++j) {
      auto d = x4c.data()[i * 4 + j];
      cout << (d > 9 ? "" : " ") << d << ' ';
    }
    cout << endl;
  }

  for (size_t i = 0; i < 12; ++i) {
    ASSERT_EQ(raw_asif_rowmajor(x4r, i), raw_asif_rowmajor(x4c, i));
    ASSERT_EQ(raw_asif_rowmajor(x4r, i), raw_asif_rowmajor(x3r, i));
    ASSERT_EQ(raw_asif_rowmajor(x4r, i), raw_asif_rowmajor(x3c, i));
    ASSERT_EQ(raw_asif_rowmajor(x4c, i), raw_asif_rowmajor(x3r, i));
    ASSERT_EQ(raw_asif_rowmajor(x4c, i), raw_asif_rowmajor(x3c, i));
    ASSERT_EQ(raw_asif_rowmajor(x3r, i), raw_asif_rowmajor(x3c, i));
  }
}
}
}