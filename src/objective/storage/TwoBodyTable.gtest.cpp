#include <gtest/gtest.h>

#include "objective/storage/TwoBodyTable.hpp"

namespace scheme {
namespace objective {
namespace storage {
namespace ritest {

using std::cout;
using std::endl;

TEST(TwoBodyTable, test) {
  TwoBodyTable<float> twob(2, 3);

  EXPECT_EQ(twob.onebody_.shape()[0], 2);
  EXPECT_EQ(twob.onebody_.shape()[1], 3);

  twob.set_onebody(0, 0, -1.0);
  twob.set_onebody(0, 1, 1.0);
  twob.set_onebody(0, 2, 2.0);
  twob.set_onebody(1, 0, 1.0);
  twob.set_onebody(1, 1, -1.0);
  twob.set_onebody(1, 2, -2.0);

  twob.init_onebody_filter(0.0);

  EXPECT_EQ(twob.all2sel_[0][0], 0);
  EXPECT_EQ(twob.all2sel_[0][1], -1);
  EXPECT_EQ(twob.all2sel_[0][2], -1);
  EXPECT_EQ(twob.all2sel_[1][0], -1);
  EXPECT_EQ(twob.all2sel_[1][1], 0);
  EXPECT_EQ(twob.all2sel_[1][2], 1);

  EXPECT_EQ(twob.nsel_[0], 1);
  EXPECT_EQ(twob.nsel_[1], 2);
  EXPECT_EQ(twob.sel2all_[0][0], 0);
  EXPECT_EQ(twob.sel2all_[0][1], -1);
  EXPECT_EQ(twob.sel2all_[0][2], -1);
  EXPECT_EQ(twob.sel2all_[1][0], 1);
  EXPECT_EQ(twob.sel2all_[1][1], 2);
  EXPECT_EQ(twob.sel2all_[1][2], -1);
}
}
}
}
}
