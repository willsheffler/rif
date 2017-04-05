#include "actor/Atom.hpp"
#include <gtest/gtest.h>
#include "eigen_types.hpp"

namespace rif {
namespace actor {
namespace test_atom {
TEST(Atom, Atom_test) {
  ASSERT_EQ(sizeof(Atom<V3f>), 16);
  Atom<V3f> a(V3f(1, 2, 3));
  ASSERT_EQ(a.position()[0], 1.0f);
  ASSERT_EQ(a.position()[1], 2.0f);
  ASSERT_EQ(a.position()[2], 3.0f);

  X3f x = X3f::Identity();
  x.translation() = V3f(9, 5, 3);

  Eigen::Map<V3f> test((float*)&a);

  Atom<V3f> a2 = Atom<V3f>(a, x);
  ASSERT_EQ(a2.position()[0], 10.0f);
  ASSERT_EQ(a2.position()[1], 7.0f);
  ASSERT_EQ(a2.position()[2], 6.0f);
}
}
}
}