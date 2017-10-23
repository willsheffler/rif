#include <gtest/gtest.h>

#include <Eigen/Geometry>
#include "geom/rand_geom.hpp"
#include "index/XformMap.hpp"

#include <random>
#include <sparsehash/dense_hash_set>
#include "util/Timer.hpp"

namespace rif {
namespace hash {
namespace xhtest {

using std::cout;
using std::endl;

typedef Eigen::Transform<double, 3, Eigen::AffineCompact> Xform;
// typedef Eigen::Affine3d Xform;

template <template <class X> class XformHash>
int get_num_ori_cells(int ori_nside, double &xcov) {
  std::mt19937 rng((unsigned int)time(0) + 7693487);
  XformHash<Xform> xh(1.0, ori_nside, 512.0);
  int n_ori_bins;
  {
    google::dense_hash_set<size_t> idx_seen;
    idx_seen.set_empty_key(std::numeric_limits<uint64_t>::max());

    int NSAMP = std::max(1000000, 500 * ori_nside * ori_nside * ori_nside);
    Xform x;
    for (int i = 0; i < NSAMP; ++i) {
      geom::rand_xform(rng, x, 512.0);
      x.translation()[0] = x.translation()[1] = x.translation()[2] = 0;
      idx_seen.insert(xh.get_key(x));
    }
    n_ori_bins = (int)idx_seen.size();
    xcov = (double)NSAMP / n_ori_bins;
  }
  return n_ori_bins;
}

TEST(XformHash_Quat_BCC7_Zorder, num_ori_cells) {
  {
    XformHash_Quat_BCC7_Zorder<Xform> xh(1.0, 2, 512.0);
    ASSERT_EQ(xh.approx_nori(), 53);
  }
  {
    XformHash_Quat_BCC7_Zorder<Xform> xh(1.0, 3, 512.0);
    ASSERT_EQ(xh.approx_nori(), 53);
  }
  {
    XformHash_Quat_BCC7_Zorder<Xform> xh(1.0, 4, 512.0);
    ASSERT_EQ(xh.approx_nori(), 181);
  }
  {
    XformHash_Quat_BCC7_Zorder<Xform> xh(1.0, 5, 512.0);
    ASSERT_EQ(xh.approx_nori(), 321);
  }
  {
    XformHash_Quat_BCC7_Zorder<Xform> xh(1.0, 6, 512.0);
    ASSERT_EQ(xh.approx_nori(), 665);
  }
  {
    XformHash_Quat_BCC7_Zorder<Xform> xh(1.0, 7, 512.0);
    ASSERT_EQ(xh.approx_nori(), 874);
  }
  {
    XformHash_Quat_BCC7_Zorder<Xform> xh(1.0, 8, 512.0);
    ASSERT_EQ(xh.approx_nori(), 1642);
  }
  {
    XformHash_Quat_BCC7_Zorder<Xform> xh(1.0, 9, 512.0);
    ASSERT_EQ(xh.approx_nori(), 1997);
  }
  {
    XformHash_Quat_BCC7_Zorder<Xform> xh(1.0, 10, 512.0);
    ASSERT_EQ(xh.approx_nori(), 2424);
  }
  {
    XformHash_Quat_BCC7_Zorder<Xform> xh(1.0, 11, 512.0);
    ASSERT_EQ(xh.approx_nori(), 3337);
  }
  {
    XformHash_Quat_BCC7_Zorder<Xform> xh(1.0, 12, 512.0);
    ASSERT_EQ(xh.approx_nori(), 4504);
  }
  {
    XformHash_Quat_BCC7_Zorder<Xform> xh(1.0, 13, 512.0);
    ASSERT_EQ(xh.approx_nori(), 5269);
  }
  {
    XformHash_Quat_BCC7_Zorder<Xform> xh(1.0, 14, 512.0);
    ASSERT_EQ(xh.approx_nori(), 6592);
  }
  {
    XformHash_Quat_BCC7_Zorder<Xform> xh(1.0, 15, 512.0);
    ASSERT_EQ(xh.approx_nori(), 8230);
  }
  {
    XformHash_Quat_BCC7_Zorder<Xform> xh(1.0, 16, 512.0);
    ASSERT_EQ(xh.approx_nori(), 10193);
  }
  {
    XformHash_Quat_BCC7_Zorder<Xform> xh(1.0, 17, 512.0);
    ASSERT_EQ(xh.approx_nori(), 11420);
  }
  {
    XformHash_Quat_BCC7_Zorder<Xform> xh(1.0, 18, 512.0);
    ASSERT_EQ(xh.approx_nori(), 14068);
  }
  {
    XformHash_Quat_BCC7_Zorder<Xform> xh(1.0, 19, 512.0);
    ASSERT_EQ(xh.approx_nori(), 16117);
  }
  {
    XformHash_Quat_BCC7_Zorder<Xform> xh(1.0, 20, 512.0);
    ASSERT_EQ(xh.approx_nori(), 19001);
  }
  {
    XformHash_Quat_BCC7_Zorder<Xform> xh(1.0, 21, 512.0);
    ASSERT_EQ(xh.approx_nori(), 21362);
  }
  {
    XformHash_Quat_BCC7_Zorder<Xform> xh(1.0, 22, 512.0);
    ASSERT_EQ(xh.approx_nori(), 25401);
  }
  {
    XformHash_Quat_BCC7_Zorder<Xform> xh(1.0, 23, 512.0);
    ASSERT_EQ(xh.approx_nori(), 29191);
  }
  {
    XformHash_Quat_BCC7_Zorder<Xform> xh(1.0, 24, 512.0);
    ASSERT_EQ(xh.approx_nori(), 33227);
  }
  {
    XformHash_Quat_BCC7_Zorder<Xform> xh(1.0, 25, 512.0);
    ASSERT_EQ(xh.approx_nori(), 37210);
  }
  {
    XformHash_Quat_BCC7_Zorder<Xform> xh(1.0, 26, 512.0);
    ASSERT_EQ(xh.approx_nori(), 41454);
  }
  {
    XformHash_Quat_BCC7_Zorder<Xform> xh(1.0, 27, 512.0);
    ASSERT_EQ(xh.approx_nori(), 45779);
  }
  {
    XformHash_Quat_BCC7_Zorder<Xform> xh(1.0, 28, 512.0);
    ASSERT_EQ(xh.approx_nori(), 51303);
  }
  {
    XformHash_Quat_BCC7_Zorder<Xform> xh(1.0, 29, 512.0);
    ASSERT_EQ(xh.approx_nori(), 57248);
  }
  {
    XformHash_Quat_BCC7_Zorder<Xform> xh(1.0, 30, 512.0);
    ASSERT_EQ(xh.approx_nori(), 62639);
  }
  {
    XformHash_Quat_BCC7_Zorder<Xform> xh(1.0, 31, 512.0);
    ASSERT_EQ(xh.approx_nori(), 69417);
  }
  {
    XformHash_Quat_BCC7_Zorder<Xform> xh(1.0, 32, 512.0);
    ASSERT_EQ(xh.approx_nori(), 76572);
  }
  {
    XformHash_Quat_BCC7_Zorder<Xform> xh(1.0, 33, 512.0);
    ASSERT_EQ(xh.approx_nori(), 83178);
  }
  {
    XformHash_Quat_BCC7_Zorder<Xform> xh(1.0, 34, 512.0);
    ASSERT_EQ(xh.approx_nori(), 92177);
  }
  {
    XformHash_Quat_BCC7_Zorder<Xform> xh(1.0, 35, 512.0);
    ASSERT_EQ(xh.approx_nori(), 99551);
  }
  {
    XformHash_Quat_BCC7_Zorder<Xform> xh(1.0, 36, 512.0);
    ASSERT_EQ(xh.approx_nori(), 108790);
  }
  {
    XformHash_Quat_BCC7_Zorder<Xform> xh(1.0, 37, 512.0);
    ASSERT_EQ(xh.approx_nori(), 117666);
  }
  {
    XformHash_Quat_BCC7_Zorder<Xform> xh(1.0, 38, 512.0);
    ASSERT_EQ(xh.approx_nori(), 127850);
  }
  {
    XformHash_Quat_BCC7_Zorder<Xform> xh(1.0, 39, 512.0);
    ASSERT_EQ(xh.approx_nori(), 138032);
  }
  {
    XformHash_Quat_BCC7_Zorder<Xform> xh(1.0, 40, 512.0);
    ASSERT_EQ(xh.approx_nori(), 149535);
  }
  {
    XformHash_Quat_BCC7_Zorder<Xform> xh(1.0, 41, 512.0);
    ASSERT_EQ(xh.approx_nori(), 159922);
  }
  {
    XformHash_Quat_BCC7_Zorder<Xform> xh(1.0, 42, 512.0);
    ASSERT_EQ(xh.approx_nori(), 171989);
  }
  {
    XformHash_Quat_BCC7_Zorder<Xform> xh(1.0, 43, 512.0);
    ASSERT_EQ(xh.approx_nori(), 183625);
  }
  {
    XformHash_Quat_BCC7_Zorder<Xform> xh(1.0, 44, 512.0);
    ASSERT_EQ(xh.approx_nori(), 196557);
  }
  {
    XformHash_Quat_BCC7_Zorder<Xform> xh(1.0, 45, 512.0);
    ASSERT_EQ(xh.approx_nori(), 209596);
  }
  {
    XformHash_Quat_BCC7_Zorder<Xform> xh(1.0, 46, 512.0);
    ASSERT_EQ(xh.approx_nori(), 226672);
  }
  {
    XformHash_Quat_BCC7_Zorder<Xform> xh(1.0, 47, 512.0);
    ASSERT_EQ(xh.approx_nori(), 239034);
  }
  {
    XformHash_Quat_BCC7_Zorder<Xform> xh(1.0, 48, 512.0);
    ASSERT_EQ(xh.approx_nori(), 253897);
  }
  {
    XformHash_Quat_BCC7_Zorder<Xform> xh(1.0, 49, 512.0);
    ASSERT_EQ(xh.approx_nori(), 271773);
  }
  {
    XformHash_Quat_BCC7_Zorder<Xform> xh(1.0, 50, 512.0);
    ASSERT_EQ(xh.approx_nori(), 288344);
  }
  {
    XformHash_Quat_BCC7_Zorder<Xform> xh(1.0, 51, 512.0);
    ASSERT_EQ(xh.approx_nori(), 306917);
  }
  {
    XformHash_Quat_BCC7_Zorder<Xform> xh(1.0, 52, 512.0);
    ASSERT_EQ(xh.approx_nori(), 324284);
  }
  {
    XformHash_Quat_BCC7_Zorder<Xform> xh(1.0, 53, 512.0);
    ASSERT_EQ(xh.approx_nori(), 342088);
  }
  {
    XformHash_Quat_BCC7_Zorder<Xform> xh(1.0, 54, 512.0);
    ASSERT_EQ(xh.approx_nori(), 364686);
  }
  {
    XformHash_Quat_BCC7_Zorder<Xform> xh(1.0, 55, 512.0);
    ASSERT_EQ(xh.approx_nori(), 381262);
  }
  {
    XformHash_Quat_BCC7_Zorder<Xform> xh(1.0, 56, 512.0);
    ASSERT_EQ(xh.approx_nori(), 405730);
  }
  {
    XformHash_Quat_BCC7_Zorder<Xform> xh(1.0, 57, 512.0);
    ASSERT_EQ(xh.approx_nori(), 427540);
  }
  {
    XformHash_Quat_BCC7_Zorder<Xform> xh(1.0, 58, 512.0);
    ASSERT_EQ(xh.approx_nori(), 450284);
  }
  {
    XformHash_Quat_BCC7_Zorder<Xform> xh(1.0, 59, 512.0);
    ASSERT_EQ(xh.approx_nori(), 472265);
  }
  {
    XformHash_Quat_BCC7_Zorder<Xform> xh(1.0, 60, 512.0);
    ASSERT_EQ(xh.approx_nori(), 498028);
  }
  {
    XformHash_Quat_BCC7_Zorder<Xform> xh(1.0, 61, 512.0);
    ASSERT_EQ(xh.approx_nori(), 521872);
  }
  {
    XformHash_Quat_BCC7_Zorder<Xform> xh(1.0, 62, 512.0);
    ASSERT_EQ(xh.approx_nori(), 547463);
  }

  // for(int i = 2; i < 64; ++i){
  //  double xcov;
  //  int num_ori_cells = get_num_ori_cells<XformHash_Quat_BCC7_Zorder>( i,
  // xcov );
  //  printf("  Quat_BCC7_Zorder %2d %7d %8.3f\n",i, num_ori_cells,
  // xcov );
  //  std::cout.flush();
  // }
}

// TEST( XformHash_bt24_BCC6, num_ori_cells ){
//  for(int i = 19; i < 65; ++i){
//    double xcov;
//    int num_ori_cells = get_num_ori_cells<XformHash_bt24_BCC6>( i,
// xcov
// );
//    printf("  bt24_BCC6 %2d %7d %8.3f\n",i, num_ori_cells,
// xcov
// );
//    std::cout.flush();
//  }
// }

template <template <class X> class XformHash>
void test_xform_hash_perf(double cart_resl, double ang_resl,
                          int const N2 = 100 * 1000, unsigned int seed = 0) {
  std::mt19937 rng((unsigned int)time(0) + seed);

  double time_key = 0.0, time_cen = 0.0;
  double cart_resl2 = cart_resl * cart_resl;
  double ang_resl2 = ang_resl * ang_resl;

  XformHash<Xform> xh(cart_resl, ang_resl, 512.0);

  std::vector<Xform> samples(N2), centers(N2);

  for (int i = 0; i < N2; ++i) geom::rand_xform(rng, samples[i], 512.0);

  util::Timer tk;
  std::vector<uint64_t> keys(N2);
  for (int i = 0; i < N2; ++i) {
    keys[i] = xh.get_key(samples[i]);
    // centers[i] = xh.get_center( keys[i] );
    // cout << endl;
  }
  time_key += (double)tk.elapsed_nano();

  util::Timer tc;
  for (int i = 0; i < N2; ++i) centers[i] = xh.get_center(keys[i]);
  time_cen += (double)tc.elapsed_nano();

  google::dense_hash_set<size_t> idx_seen;
  idx_seen.set_empty_key(std::numeric_limits<uint64_t>::max());
  for (int i = 0; i < N2; ++i) idx_seen.insert(keys[i]);

  double covrad = 0, max_dt = 0, max_da = 0;
  for (int i = 0; i < N2; ++i) {
    Xform l = centers[i].inverse() * samples[i];
    double dt = l.translation().norm();
    Eigen::Matrix3d m;
    for (int k = 0; k < 9; ++k) m.data()[k] = l.data()[k];
    // cout << m << endl;
    // cout << l.rotation() << endl;
    double da = Eigen::AngleAxisd(m).angle() * 180.0 / M_PI;
    // double da = Eigen::AngleAxisd(l.rotation()).angle()*180.0/M_PI;
    double err = sqrt(da * da / ang_resl2 * cart_resl2 + dt * dt);
    covrad = fmax(covrad, err);
    max_dt = fmax(max_dt, dt);
    max_da = fmax(max_da, da);
    ASSERT_LT(dt, cart_resl);
    ASSERT_LT(da, ang_resl);
  }
  ASSERT_GT(max_dt * 1.25, cart_resl);
  ASSERT_GT(max_da * 1.4, ang_resl);  // multiplier of 1.4??

  double tot_cell_vol = covrad * covrad * covrad * covrad * covrad * covrad *
                        xh.approx_nori() / (cart_resl * cart_resl * cart_resl);
  printf(
      " %5.3f/%5.1f cr %5.3f dt %5.3f da %6.3f x2k: %7.3fns k2x: %7.3fns "
      "%9.3f "
      "%7lu\n",
      cart_resl, ang_resl, covrad, max_dt, max_da / ang_resl, time_key / N2,
      time_cen / N2, tot_cell_vol, xh.approx_nori());

  // cout << " rate " << N1*N2/time_key << "  " << N1*N2/time_cen << endl;
}

TEST(XformHash, DISABLED_XformHash_Quat_BCC7_Zorder) {
  unsigned int s = 0;
  int N = 10 * 1000;
#ifdef SCHEME_BENCHMARK
  N = 1 * 1000 * 1000;
#endif
  cout << "  Quat_BCC7_Zorder";
  test_xform_hash_perf<XformHash_Quat_BCC7_Zorder>(4.00, 30.0, N, ++s);
  cout << "  Quat_BCC7_Zorder";
  test_xform_hash_perf<XformHash_Quat_BCC7_Zorder>(2.00, 20.0, N, ++s);
  cout << "  Quat_BCC7_Zorder";
  test_xform_hash_perf<XformHash_Quat_BCC7_Zorder>(1.00, 15.0, N, ++s);
  cout << "  Quat_BCC7_Zorder";
  test_xform_hash_perf<XformHash_Quat_BCC7_Zorder>(0.50, 10.0, N, ++s);
  cout << "  Quat_BCC7_Zorder";
  test_xform_hash_perf<XformHash_Quat_BCC7_Zorder>(0.25, 5.0, N, ++s);
  cout << "  Quat_BCC7_Zorder";
  test_xform_hash_perf<XformHash_Quat_BCC7_Zorder>(0.11, 3.3, N, ++s);
}

TEST(XformHash, XformHash_bt24_BCC6) {
  unsigned int s = 0;
  int N = 10 * 1000;
#ifdef SCHEME_BENCHMARK
  N = 1 * 1000 * 1000;
#endif
  cout << "  bt24_BCC6";
  test_xform_hash_perf<XformHash_bt24_BCC6>(4.00, 30.0, N, ++s);
  cout << "  bt24_BCC6";
  test_xform_hash_perf<XformHash_bt24_BCC6>(2.00, 20.0, N, ++s);
  cout << "  bt24_BCC6";
  test_xform_hash_perf<XformHash_bt24_BCC6>(1.00, 15.0, N, ++s);
  cout << "  bt24_BCC6";
  test_xform_hash_perf<XformHash_bt24_BCC6>(0.50, 10.0, N, ++s);
  cout << "  bt24_BCC6";
  test_xform_hash_perf<XformHash_bt24_BCC6>(0.25, 5.0, N, ++s);
  cout << "  bt24_BCC6";
  test_xform_hash_perf<XformHash_bt24_BCC6>(0.11, 3.3, N, ++s);
}

// TEST( XformHash, XformHash_bt24_BCC3 ){
//  unsigned int s = 0;
//  int N = 10*1000;
//  #ifdef SCHEME_BENCHMARK
//  N = 1*1000*1000;
//  #endif
//  cout << "  bt24_BCC3"; test_xform_hash_perf< XformHash_bt24_BCC3
// >( 4.00 , 30.0, N, ++s );
//  cout << "  bt24_BCC3"; test_xform_hash_perf< XformHash_bt24_BCC3
// >( 2.00 , 20.0, N, ++s );
//  cout << "  bt24_BCC3"; test_xform_hash_perf< XformHash_bt24_BCC3
// >( 1.00 , 15.0, N, ++s );
//  cout << "  bt24_BCC3"; test_xform_hash_perf< XformHash_bt24_BCC3
// >( 0.50 , 10.0, N, ++s );
//  cout << "  bt24_BCC3"; test_xform_hash_perf< XformHash_bt24_BCC3
// >( 0.25 ,  5.0, N, ++s );
//  cout << "  bt24_BCC3"; test_xform_hash_perf< XformHash_bt24_BCC3
// >( 0.11 ,  3.3, N, ++s );
// }

// TEST( XformHash, preformance ){

//  cout << "  Quat_BCC7        "; test_xform_hash_perf<
// XformHash_Quat_BCC7         >( 2.0 , 1 );
//  cout << "  Quat_BCC7_Zorder "; test_xform_hash_perf<
// XformHash_Quat_BCC7_Zorder  >( 2.0 , 2 );
//  cout << "  bt24_BCC3        "; test_xform_hash_perf<
// XformHash_bt24_BCC3         >( 2.0 , 3 );
//  cout << "  bt24_BCC3_Zorder "; test_xform_hash_perf<
// XformHash_bt24_BCC3_Zorder  >( 2.0 , 4 );
//  cout << "  bt24_Cubic_Zorder"; test_xform_hash_perf<
// XformHash_bt24_Cubic_Zorder >( 2.0 , 5 );

//  cout << "  Quat_BCC7        "; test_xform_hash_perf<
// XformHash_Quat_BCC7         >( 0.25 , 1 );
//  cout << "  Quat_BCC7_Zorder "; test_xform_hash_perf<
// XformHash_Quat_BCC7_Zorder  >( 0.25 , 2 );
//  cout << "  bt24_BCC3        "; test_xform_hash_perf<
// XformHash_bt24_BCC3         >( 0.25 , 3 );
//  cout << "  bt24_BCC3_Zorder "; test_xform_hash_perf<
// XformHash_bt24_BCC3_Zorder  >( 0.25 , 4 );
//  cout << "  bt24_Cubic_Zorder"; test_xform_hash_perf<
// XformHash_bt24_Cubic_Zorder >( 0.25 , 5 );
// }

TEST(XformHash, XformHash_Quat_BCC7_Zorder_cart_shift) {
  std::mt19937 rng((unsigned int)time(0) + 23908457);
  std::uniform_real_distribution<> runif;
  int NSAMP = 3;
#ifdef SCHEME_BENCHMARK
  NSAMP = 20;
#endif
  using namespace Eigen;
  for (int i = 0; i < 100; ++i) {
    XformHash_Quat_BCC7_Zorder<Xform> h(0.1 + runif(rng), 11.0,
                                        300.0 * runif(rng) + 100.0);
    Xform m = Xform::Identity();
    uint64_t mkey = h.get_key(m);
    ASSERT_EQ(m.translation(), Vector3d(0, 0, 0));
    Vector3d v;
    ASSERT_FALSE(mkey & 1);  // middle key is even
    for (int x = -NSAMP; x <= NSAMP; ++x) {
      for (int y = -NSAMP; y <= NSAMP; ++y) {
        for (int z = -NSAMP; z <= NSAMP; ++z) {
          // Odd center cell
          v = h.get_center(h.cart_shift_key(mkey, x, y, z)).translation() /
              h.cart_spacing();
          ASSERT_LE((Vector3d(x, y, z) - v).squaredNorm(), 0.00001);
          // Even center cell
          v = h.get_center(h.cart_shift_key(mkey | 1, x, y, z)).translation() /
              h.cart_spacing();
          ASSERT_LE((Vector3d(x + 0.5, y + 0.5, z + 0.5) - v).squaredNorm(),
                    0.00001);
        }
      }
    }
  }
  XformHash_Quat_BCC7_Zorder<Xform> h(0.1 + runif(rng), 11.0,
                                      300.0 * runif(rng) + 100.0);
  ASSERT_EQ(h.cart_shift_key(3450, 0, 0, 0, 0), 3450);
  ASSERT_EQ(h.cart_shift_key(3450, 0, 0, 0, 1), 3451);
  ASSERT_EQ(h.cart_shift_key(3451, 0, 0, 0, 0), 3451);
  ASSERT_EQ(h.cart_shift_key(3451, 0, 0, 0, 1), 3450);
}

TEST(XformHash, XformAngHash_bt24_BCC) {
  using XH = XformAngHash_bt24_BCC6<Xform>;
  int N1 = 5;
  int N2 = 10000;
  std::mt19937 rng((unsigned int)time(0) + 34979234);
  std::uniform_real_distribution<> runif;
  double time_key = 0.0, time_cen = 0.0;
  for (int i = 0; i < N1; ++i) {
    float phi_resl = 5.0 + runif(rng) * 20.0;
    float cart_resl = 0.1 + runif(rng);
    float ang_resl = 5.0 + runif(rng) * 20.0;
    float cart_bound = 100.0 * runif(rng) + 100.0;
    float cart_resl2 = cart_resl * cart_resl;
    float ang_resl2 = ang_resl * ang_resl;
    XH xh(phi_resl, cart_resl, ang_resl, cart_bound);

    std::vector<std::pair<Xform, float>> samples(N2), centers(N2);
    for (int i = 0; i < N2; ++i) {
      geom::rand_xform(rng, samples[i].first, cart_bound * 0.9);
      samples[i].second = runif(rng) * 360.0 - 180.0;
    }

    util::Timer tk;
    std::vector<uint64_t> keys(N2);
    for (int i = 0; i < N2; ++i) {
      keys[i] = xh.get_key(samples[i]);
    }
    time_key += (double)tk.elapsed_nano();

    util::Timer tc;
    for (int i = 0; i < N2; ++i) {
      centers[i] = xh.get_center(keys[i]);
      // auto k = xh.get_key(centers[i]);
      // ASSERT_EQ(k, keys[i]);
    }
    time_cen += (double)tc.elapsed_nano();

    google::dense_hash_set<size_t> idx_seen;
    idx_seen.set_empty_key(std::numeric_limits<uint64_t>::max());
    for (int i = 0; i < N2; ++i) idx_seen.insert(keys[i]);

    double covrad = 0, max_dt = 0, max_da = 0, max_dp = 0;
    for (int i = 0; i < N2; ++i) {
      Xform l = centers[i].first.inverse() * samples[i].first;
      float dp = fabs(centers[i].second - samples[i].second);
      double dt = l.translation().norm();
      double da = Eigen::AngleAxisd(l.linear()).angle() * 180.0 / M_PI;
      // double da = Eigen::AngleAxisd(l.rotation()).angle()*180.0/M_PI;
      // double err =
      // sqrt(da * da / ang_resl2 * cart_resl2 + dt * dt + dp / phi_resl);
      // covrad = fmax(covrad, err);
      max_dt = fmax(max_dt, dt);
      max_da = fmax(max_da, da);
      max_dp = fmax(max_dp, dp);
      ASSERT_LT(dt, cart_resl);
      ASSERT_LT(da, ang_resl);
      ASSERT_LT(dp, phi_resl);
    }
    ASSERT_GT(max_dt * 1.25, cart_resl);
    ASSERT_GT(max_da * 1.6, ang_resl);  // multiplier of 1.4??
    ASSERT_GT(max_dp * 1.2, phi_resl);
  }
}

TEST(XformHash, Xform2AngHash_bt24_BCC6) {
  using XH = Xform2AngHash_bt24_BCC6<Xform>;
  int N1 = 5;
  int N2 = 10000;
  std::mt19937 rng((unsigned int)time(0) + 34979234);
  std::uniform_real_distribution<> runif;
  double time_key = 0.0, time_cen = 0.0;
  for (int i = 0; i < N1; ++i) {
    float phi_resl = 5.0 + runif(rng) * 20.0;
    float cart_resl = 0.1 + runif(rng);
    float ang_resl = 5.0 + runif(rng) * 20.0;
    float cart_bound = 100.0 * runif(rng) + 100.0;
    float cart_resl2 = cart_resl * cart_resl;
    float ang_resl2 = ang_resl * ang_resl;
    XH xh(phi_resl, cart_resl, ang_resl, cart_bound);

    std::vector<std::tuple<Xform, float, float>> samples(N2), centers(N2);
    for (int i = 0; i < N2; ++i) {
      geom::rand_xform(rng, std::get<0>(samples[i]), cart_bound * 0.9);
      std::get<1>(samples[i]) = runif(rng) * 360.0 - 180.0;
      std::get<2>(samples[i]) = runif(rng) * 360.0 - 180.0;
    }

    util::Timer tk;
    std::vector<uint64_t> keys(N2);
    for (int i = 0; i < N2; ++i) {
      keys[i] = xh.get_key(samples[i]);
    }
    time_key += (double)tk.elapsed_nano();

    util::Timer tc;
    for (int i = 0; i < N2; ++i) {
      centers[i] = xh.get_center(keys[i]);
      // auto k = xh.get_key(centers[i]);
      // ASSERT_EQ(k, keys[i]);
    }
    time_cen += (double)tc.elapsed_nano();

    google::dense_hash_set<size_t> idx_seen;
    idx_seen.set_empty_key(std::numeric_limits<uint64_t>::max());
    for (int i = 0; i < N2; ++i) idx_seen.insert(keys[i]);

    double covrad = 0, max_dt = 0, max_da = 0, max_dp = 0;
    for (int i = 0; i < N2; ++i) {
      Xform l = std::get<0>(centers[i]).inverse() * std::get<0>(samples[i]);
      float dp = fabs(std::get<1>(centers[i]) - std::get<1>(samples[i]));
      float dp2 = fabs(std::get<2>(centers[i]) - std::get<2>(samples[i]));
      double dt = l.translation().norm();
      double da = Eigen::AngleAxisd(l.linear()).angle() * 180.0 / M_PI;
      // double da = Eigen::AngleAxisd(l.rotation()).angle()*180.0/M_PI;
      // double err =
      // sqrt(da * da / ang_resl2 * cart_resl2 + dt * dt + dp / phi_resl);
      // covrad = fmax(covrad, err);
      max_dt = fmax(max_dt, dt);
      max_da = fmax(max_da, da);
      max_dp = fmax(max_dp, dp);
      ASSERT_LT(dt, cart_resl);
      ASSERT_LT(da, ang_resl);
      ASSERT_LT(dp, phi_resl);
    }
    ASSERT_GT(max_dt * 1.25, cart_resl);
    ASSERT_GT(max_da * 1.6, ang_resl);  // multiplier of 1.4??
    ASSERT_GT(max_dp * 1.2, phi_resl);
  }
}

// end ns
}
}
}
