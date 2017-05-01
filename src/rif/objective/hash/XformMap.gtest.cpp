#include <gtest/gtest.h>

#include <Eigen/Geometry>
#include "geom/rand_geom.hpp"
#include "objective/hash/XformMap.hpp"

#include <sparsehash/dense_hash_set>

#include <random>
#include "util/Timer.hpp"

#include <fstream>

namespace rif {
namespace objective {
namespace hash {
namespace xmtest {

using std::cout;
using std::endl;

typedef Eigen::Transform<double, 3, Eigen::AffineCompact> Xform;
// typedef Eigen::Affine3d Xform;

TEST(XformMap, stores_correctly) {
  int NSAMP = 100000;

  std::mt19937 rng((unsigned int)time(0) + 296720384);
  std::uniform_real_distribution<> runif;

  XformMap<Xform, double> xmap(0.5, 10.0);
  std::vector<std::pair<Xform, double>> dat;
  for (int i = 0; i < NSAMP; ++i) {
    Xform x;
    geom::rand_xform(rng, x, 256.0);
    double val = runif(rng);
    xmap.insert(x, val);
    dat.push_back(std::make_pair(x, val));
  }

  XformMap<Xform, double> const &xmap_test(xmap);

  util::Timer t;

  for (int i = 0; i < dat.size(); ++i) {
    Xform const &x = dat[i].first;
    double const v = dat[i].second;
    EXPECT_EQ(xmap_test[x], v);
    // cout << x.translation().transpose() << " " << v << endl;
  }
  cout << "XformMap " << NSAMP
       << " lookup rate: " << (double)NSAMP / t.elapsed() << " /sec ";

  // { // no way to check if stream in binary!
  //  std::cout << "following failure message is expected" << std::endl;
  //  std::ofstream out("test.sxm" );// , std::ios::binary );
  //  ASSERT_FALSE( xmap.save( out, "foo" ) );
  //  out.close();
  // }
  std::ofstream out("test.sxm", std::ios::binary);
  ASSERT_TRUE(xmap.save(out, "foo"));
  out.close();

  XformMap<Xform, double> xmap_loaded;
  std::ifstream in("test.sxm", std::ios::binary);
  ASSERT_TRUE(xmap_loaded.load(in));
  in.close();

  ASSERT_EQ(xmap.cart_resl_, xmap_loaded.cart_resl_);
  ASSERT_EQ(xmap.ang_resl_, xmap_loaded.ang_resl_);
  for (int i = 0; i < dat.size(); ++i) {
    Xform const &x = dat[i].first;
    double const v = dat[i].second;
    ASSERT_EQ(xmap.hasher_.get_key(x), xmap_loaded.hasher_.get_key(x));
    ASSERT_EQ(xmap_loaded[x], v);
    // cout << x.translation().transpose() << " " << v << endl;
  }
}

double get_ident_lever_dis(Xform x, double lever_dis) {
  util::SimpleArray<7, double> x_lever_coord;
  x_lever_coord[0] = x.translation()[0];
  x_lever_coord[1] = x.translation()[1];
  x_lever_coord[2] = x.translation()[2];
  Eigen::Matrix<double, 3, 3> rot;
  get_transform_rotation(x, rot);
  Eigen::Quaternion<double> q(rot);
  bool neg = q.w() < 0.0;
  x_lever_coord[3] = (neg ? -q.w() : q.w()) * 2.0 * lever_dis;
  x_lever_coord[4] = (neg ? -q.x() : q.x()) * 2.0 * lever_dis;
  x_lever_coord[5] = (neg ? -q.y() : q.y()) * 2.0 * lever_dis;
  x_lever_coord[6] = (neg ? -q.z() : q.z()) * 2.0 * lever_dis;
  util::SimpleArray<7, double> ident(0, 0, 0, lever_dis * 2.0, 0, 0, 0);
  return (x_lever_coord - ident).norm();
}

TEST(XformMap, insert_sphere) {
  int NSAMP2 = 10000;

  // typedef XformMap< Xform, double, 0 > XMap;
  typedef XformMap<Xform, double> XMap;
  std::mt19937 rng((unsigned int)time(0) + 3457820);
  std::uniform_real_distribution<> runif;
  Xform x;
  double cart_resl = 1.0;
  double lever = 3.0;
  double ang_resl = cart_resl / lever * 180.0 / M_PI;
  XMap xmap(1.00, ang_resl);
  double rad = 3.0;
  cout << "cart_resl " << cart_resl << " ang_resl " << ang_resl << " lever "
       << lever << " sphere rad " << rad << endl;
  double angrad = rad / lever * 180.0 / M_PI;
  double quatrad = numeric::deg2quat(angrad);
  geom::rand_xform(rng, x, 256.0);
  XformHashNeighbors<XMap::Hasher> nbcache(rad, angrad, xmap.hasher_, 500.0);
  int nbitercount = xmap.insert_sphere(x, rad, lever, 12345.0, nbcache);
  cout << nbitercount << " " << xmap.count(12345.0) << " "
       << xmap.size() - (float)xmap.count(0) << " "
       << (float)xmap.count(0) / xmap.size() << " " << xmap.map_.size() << endl;
  int n_cart_fail = 0, n_rot_fail = 0, n_both_fail = 0;
  int n_lever_false_pos = 0, n_lever_false_neg = 0, n_within = 0, n_without = 0;
  for (int i = 0; i < NSAMP2; ++i) {
    Xform p;
    geom::rand_xform_quat(rng, p, rad, 0.0);
    if (xmap[x * p] != 12345.0) ++n_cart_fail;

    geom::rand_xform_quat(rng, p, 0.0, quatrad);
    if (xmap[x * p] != 12345.0) ++n_rot_fail;

    double split1 = runif(rng);
    double split2 = sqrt(1.0 - split1 * split1);
    geom::rand_xform_quat(rng, p, split1 * rad, split2 * quatrad);
    if (xmap[x * p] != 12345.0) ++n_both_fail;

    geom::rand_xform_quat(rng, p, 1.5 * split1 * rad, 1.5 * split2 * quatrad);
    // geom::rand_xform_quat(rng, p, rad, quatrad );
    if (get_ident_lever_dis(p, lever) < rad) {
      ++n_within;
      if (xmap[x * p] != 12345.0) ++n_lever_false_neg;
    } else {
      ++n_without;
      if (xmap[x * p] == 12345.0) ++n_lever_false_pos;
    }
  }
  cout << "CART FAIL FRAC " << (float)n_cart_fail / NSAMP2 << endl;
  cout << "ROT  FAIL FRAC " << (float)n_rot_fail / NSAMP2 << endl;
  cout << "BOTH FAIL FRAC " << (float)n_both_fail / NSAMP2 << endl;

  cout << "FALSE POS " << (float)n_lever_false_pos / n_without << endl;
  cout << "FALSE NEG " << (float)n_lever_false_neg / n_within << ",  FRAC "
       << (float)n_within / NSAMP2 << endl;

  ASSERT_LT((float)n_cart_fail / NSAMP2, 0.03);
  ASSERT_LT((float)n_rot_fail / NSAMP2, 0.05);
  ASSERT_LT((float)n_both_fail / NSAMP2, 0.02);
  ASSERT_LT((float)n_lever_false_pos / n_without, 0.30);
  ASSERT_LT((float)n_lever_false_neg / n_within, 0.015);
}

TEST(XformMap, test_bt24_bcc6) {
  typedef Eigen::Transform<double, 3, Eigen::AffineCompact> EigenXform;
  typedef rif::objective::hash::XformMap<EigenXform, double,
                                         XformHash_bt24_BCC6>
      XMap;

  XMap xmap(1.0, 10.0);
}

TEST(XformMap, DISABLED_test_float_double) {
  typedef Eigen::Transform<double, 3, Eigen::AffineCompact> EigenXformD;
  typedef rif::objective::hash::XformMap<EigenXformD, double,
                                         XformHash_bt24_BCC6>
      XMapD;
  typedef Eigen::Transform<float, 3, Eigen::AffineCompact> EigenXformF;
  typedef rif::objective::hash::XformMap<EigenXformF, double,
                                         XformHash_bt24_BCC6>
      XMapF;

  ASSERT_TRUE(false);
}
}
}
}
}
