#include <gtest/gtest.h>

#include "rif/geom/rand_geom.hpp"
#include "rif/index/stripe_index_3d.hpp"
#include "rif/util/Timer.hpp"

namespace rif {
namespace index {
namespace test_xyzhash {

using rif::index::stripe_index_3d;
using namespace std;
using F = float;

TEST(stripe_index_3d, test_nopayload_with_timing) {
  int const Niter = 10;
  int const Ntest = 100;
  int const Npts = 1000;
  float radius = 0.1;  // should be avg nbr ~3?
  int64_t num0 = 0, numtot = 0, numnbr = 0;
  double btime = 0, htime = 0, ctime = 0;
  for (int iter = 0; iter < Niter; ++iter) {
    if (iter > 0 && iter % 1000 == 0)
      cout << iter << " " << numtot << " " << num0 << endl;
    auto pts = geom::rand_box_n<F>(Npts);
    stripe_index_3d<V3f> h(radius, pts);
    h.sanity_check();
    vector<int> hcount(Ntest), bcount(Ntest), contact(Ntest);
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
    // cout << "hcount " << hcount << " bcount " << bcount << endl;
    for (int itest = 0; itest < Ntest; ++itest) {
      ASSERT_EQ(hcount[itest], bcount[itest]);
      ASSERT_EQ(contact[itest], bcount[itest] != 0);
      num0 += hcount[itest] == 0;
      numtot++;
      numnbr += hcount[itest];
    }
  }
  float avg_nbr = float(numnbr) / float(numtot);
  cout << "             " << numtot << " " << num0
       << " nbcount speedup: " << btime / htime << "x (" << btime / ctime
       << "x)"
       << " avgNnbrs: " << avg_nbr
       << " Npts/(avgNnbrs): " << float(Npts) / avg_nbr << endl;
}

TEST(stripe_index_3d, test_with_payload) {
  int const Niter = 100;
  int const Ntest = 100;
  int const Npts = 100;
  float radius = 0.25;
  float fudge_factor = 1.0;  // 1.0 == no fudge!
  double nfail = 0, ntot = 0;
  for (int iter = 0; iter < Niter; ++iter) {
    auto pts = geom::rand_box_n<F>(Npts);
    for (auto& v : pts) v = v + V3<F>(1 / sqrt(3), 1 / sqrt(3), 1 / sqrt(3));
    // pts[0] = V3f(0, 0, 0);
    stripe_index_3d<V3f, V3f> h(radius, pts, pts);
    h.sanity_check();
    auto pts_test = geom::rand_box_n<F>(Ntest);
    for (auto& v : pts_test) v = 1.1 * v - V3<F>(0.05, 0.05, 0.05);
    for (int itest = 0; itest < Ntest; ++itest) {
      auto q = pts_test[itest];
      int nbcount = h.nbcount(q);
      ASSERT_EQ(nbcount, h.brute_nbcount(q));

      vector<V3f> neighboring_points = h.neighboring_points(q);
      ASSERT_EQ(nbcount, neighboring_points.size());
      for (V3f n : neighboring_points) {
        float dist = (n - q).norm();
        if (dist > radius) ++nfail;
        ++ntot;
        ASSERT_LE(dist, radius * fudge_factor);
      }

      vector<V3f> neighboring_payloads = h.neighboring_payloads(q);
      ASSERT_EQ(nbcount, neighboring_payloads.size());
      for (V3f n : neighboring_payloads) {
        ASSERT_LE((n - q).norm(), radius * fudge_factor);
      }

      vector<pair<V3f, V3f>> neighboring_values = h.neighboring_values(q);
      ASSERT_EQ(nbcount, neighboring_values.size());
      for (auto n : neighboring_values) {
        ASSERT_LE((n.first - q).norm(), radius * fudge_factor);
        ASSERT_LE((n.second - q).norm(), radius * fudge_factor);
        ASSERT_LE((n.first - n.second).norm(), 0.0001);
      }
    }
    if (nfail) cout << nfail / ntot << ' ' << nfail << ' ' << ntot << endl;
  }
}
}
}
}
