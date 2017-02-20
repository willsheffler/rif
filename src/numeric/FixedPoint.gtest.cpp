#include <gtest/gtest.h>

#include "numeric/FixedPoint.hpp"

namespace scheme {
namespace numeric {
namespace test {

using std::cout;
using std::endl;

TEST(FixedPoint, test) {
  typedef FixedPoint<17> F17;
  typedef FixedPoint<-17> Fm17;
  Fm17 f;

  for (float d = 0; d > -15.0; d -= 0.1) {
    f = d;
    ASSERT_FLOAT_EQ((float)f, (int)(d * -17.0) / -17.0);
  }

  f = -999.0;
  ASSERT_FLOAT_EQ((float)f, -255.0 / 17.0);

  f = 1.0;
  ASSERT_FLOAT_EQ((float)f, 0.0);
}
}
}
}
