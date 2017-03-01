#include <gtest/gtest.h>

#include "search/HackPack.hpp"

namespace rif {
namespace search {
namespace hptest {

using std::cout;
using std::endl;

typedef float Float;

TEST(HackPack, create_empty_packer) {
  ::rif::objective::storage::TwoBodyTable<float> twob(1, 1);
  HackPackOpts opts;

  HackPack packer(twob, opts, 0);
}
}
}
}
