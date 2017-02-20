#include <gtest/gtest.h>

#include "nest/NEST.hpp"
#include "nest/pmap/OriTransMap.hpp"

namespace scheme {
namespace nest {
namespace pmap {

using std::cout;
using std::endl;

TEST(OriTransMap, basic_test) {
  typedef Eigen::Transform<double, 3, Eigen::AffineCompact> EigenXform;

  typedef NEST<6, EigenXform, OriTransMap> Nest;

  Nest nest(999.0, 0.0, 100.0, 1);

  std::cout << "THIS IS NO TEST!" << std::endl;
  for (int i = 0; i < nest.size(1) / 192; ++i) {
    bool valid = nest.set_state(i, 1);
    if (!valid) continue;
    EigenXform x = nest.value();
    Eigen::AngleAxisd aa(x.rotation());
    std::cout << i << " " << x.translation().transpose() << " " << aa.angle()
              << " " << aa.axis().transpose() << std::endl;
  }
}

TEST(OriTransMap, lookups_multicell_trans) {
  typedef Eigen::Transform<double, 3, Eigen::AffineCompact> EigenXform;

  typedef NEST<6, EigenXform, OriTransMap> Nest;

  Nest nest(80.0, -512.0, 512.0, 1024);
  // cout << nest.size(0) / (1024ull*1024*1024) << endl;

  int const resl = 2;
  int end = std::min((int)1000000, (int)nest.size(resl));
  for (int i = 0; i < end; ++i) {
    EigenXform x;
    if (!nest.get_state(i, resl, x)) continue;
    ASSERT_EQ(nest.get_index(x, resl), i);
  }
}

TEST(OriTransMap, lookups_multicell_ori) {
  typedef Eigen::Transform<double, 3, Eigen::AffineCompact> EigenXform;

  typedef NEST<6, EigenXform, OriTransMap> Nest;

  Nest nest(10.0, 0.0, 1.0, 1);

  int const resl = 0;
  int end = std::min((int)1000000, (int)nest.size(resl));
  int nfail = 0;
  for (int i = 0; i < end; ++i) {
    EigenXform x;
    if (!nest.get_state(i, resl, x)) continue;
    nfail += (nest.get_index(x, resl) != i);
  }

  ASSERT_LE(nfail * 1. / end, 0.1);
}

TEST(OriTransMap, name) {
  typedef Eigen::Transform<double, 3, Eigen::AffineCompact> EigenXform;

  typedef NEST<6, EigenXform, OriTransMap> Nest;

  ASSERT_EQ(

      Nest::pmap_name(),

      "OriTransMap< TetracontoctachoronMap, ScaleMap<3> >"

      );
}
}
}
}
