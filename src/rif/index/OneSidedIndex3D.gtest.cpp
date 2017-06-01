#include <gtest/gtest.h>

#include "rif/geom/rand_geom.hpp"
#include "rif/index/OneSidedIndex3D.hpp"
#include "rif/util/Timer.hpp"

namespace rif {
namespace index {
namespace test_xyzhash {

using rif::index::xyzStripeHash;

using F = float;

TEST(OneSidedIndex3D, test_vs_brute_force) {
  int const Niter = 10;
  int const Ntest = 1000;
  int const Npts = 1000;
  int64_t num0 = 0, numtot = 0;
  double btime = 0, htime = 0, ctime = 0;
  for (int iter = 0; iter < Niter; ++iter) {
    if (iter > 0 && iter % 1000 == 0)
      std::cout << iter << " " << numtot << " " << num0 << std::endl;
    auto pts = geom::rand_box_n<F>(Npts);
    xyzStripeHash<V3f> h(0.07, pts);
    // h.sanity_check();
    std::vector<int> hcount(Ntest), bcount(Ntest), contact(Ntest);
    auto pts_test = geom::rand_box_n<F>(Ntest);
    for (auto& v : pts_test) v = 1.1 * v - V3<F>(0.05, 0.05, 0.05);
    util::Timer hcount_timer;
    for (int itest = 0; itest < Ntest; ++itest)
      hcount[itest] = h.nbcount(pts_test[itest]);
    hcount_timer.stop();
    util::Timer bcount_timer;
    for (int itest = 0; itest < Ntest; ++itest)
      bcount[itest] = h.brute_nbcount(pts_test[itest]);
    bcount_timer.stop();
    util::Timer contact_timer;
    for (int itest = 0; itest < Ntest; ++itest)
      contact[itest] = h.contact(pts_test[itest]);
    contact_timer.stop();
    htime += hcount_timer.elapsed();
    btime += bcount_timer.elapsed();
    ctime += contact_timer.elapsed();
    // std::cout << "hcount " << hcount << " bcount " << bcount << std::endl;
    for (int itest = 0; itest < Ntest; ++itest) {
      ASSERT_EQ(hcount[itest], bcount[itest]);
      ASSERT_EQ(contact[itest], bcount[itest] != 0);
      num0 += hcount[itest] == 0;
      numtot++;
    }
  }
  std::cout << "             " << numtot << " " << num0
            << " nbcount speedup: " << btime / htime << "x (" << btime / ctime
            << "x)" << std::endl;
}
}
}
}
