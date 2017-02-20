#include <gtest/gtest.h>

#include "nest/NEST.hpp"
#include "nest/pmap/QuaternionMap.hpp"

#include <boost/format.hpp>
#include <boost/lexical_cast.hpp>
#include <random>

#include "io/dump_pdb_atom.hpp"

#include <Eigen/Geometry>
#include <fstream>

namespace scheme {
namespace nest {
namespace pmap {

using std::cout;
using std::endl;

TEST(QuaternionMap, test_cell_validity_check) {
  ASSERT_EQ((int)1.5, 1);  // check round towards 0 ... IEEE?
  ASSERT_EQ((int)-1.5, -1);
  typedef QuaternionMap<4> MapType;
  typedef MapType::Params P;
  BOOST_STATIC_ASSERT(
      (util::meta::has_subscript_oper<P, double &, size_t>::value));
  BOOST_STATIC_ASSERT(
      (util::meta::has_const_subscript_oper<P, double const &, size_t>::value));

  MapType qmap;
  MapType::ValueType val;
  P q;

  // q = P( 0.00, 0.00, 0.00, 0.00 ); ASSERT_TRUE ( qmap.params_to_value(
  // q/2.0+0.5, 0, 2, val ) );
  // q = P(-0.01, 0.00, 0.00, 0.00 ); ASSERT_FALSE( qmap.params_to_value(
  // q/2.0+0.5, 0, 2, val ) );
  // q = P( 0.01,-0.01, 0.01,-0.01 ); ASSERT_TRUE ( qmap.params_to_value(
  // q/2.0+0.5, 0, 2, val ) );
  // q = P( 0.01,-0.01, 0.01,-0.01 ); ASSERT_FALSE( qmap.params_to_value(
  // q/2.0+0.5, 0, 3, val ) );

  NEST<4, P, QuaternionMap> nest;

  // test that a reasonable fraction of values are invalid
  int count = 0;
  for (size_t resl = 2; resl < 6; ++resl) {
    int nvalid = 0;
    for (size_t i = 0; i < nest.size(resl); ++i) {
      if (nest.set_state(i, resl)) {
        ++nvalid;
      }
      ++count;
    }
    double empirical_frac_valid = (double)nvalid / (double)nest.size(resl);
    // pi*pi is hyper area of 1/2 3sphere
    // 4.0*1<<resl is 2 * cell width
    // 16.0 is vol of 4-cube
    double analytic_frac_hypervolume_aprox =
        3.14159 * 3.14159 * (4.0 / (1 << resl)) / 16.0;
    cout << "resl " << resl << " frac valid " << empirical_frac_valid << " vs. "
         << analytic_frac_hypervolume_aprox
         << " approx cell-width hyper-shell volume fraction" << endl;
    ASSERT_LT(empirical_frac_valid, analytic_frac_hypervolume_aprox);
    if (resl > 2)
      ASSERT_GT(empirical_frac_valid, analytic_frac_hypervolume_aprox / 2.0);
  }
  cout << "Num QuaternionMap checks: " << (double)count / 1000000.0 << endl;

  std::mt19937 rng((unsigned int)time(0));
  std::normal_distribution<> gauss;
  std::uniform_real_distribution<> uniform;

  // test that the indices of all valid quats are valid for set_state
  for (int i = 0; i < 10000; ++i) {
    q = P(fabs(gauss(rng)), gauss(rng), gauss(rng), gauss(rng));
    q = q / q.norm();
    // cout << q << endl;
    for (size_t resl = 0; resl < 15; ++resl) {
      size_t i_should_be_valid = nest.get_index(q, resl);
      ASSERT_TRUE(nest.set_state(i_should_be_valid, resl));
    }
  }
}
TEST(QuaternionMap, covering) {
  std::mt19937 rng((unsigned int)time(0));
  std::normal_distribution<> gauss;
  std::uniform_real_distribution<> uniform;

  // 	// 	1                   16 119.1666
  // 	//  2                  256 59.69261
  // 	//  3                 4096 27.35909
  // 	//  4                65536 13.58733
  // 	//  5              1048576 6.926648
  // 	//  6             16777216 3.517071
  // 	//  7            268435456 1.735485
  // 	//  8           4294967296 0.8686822
  // 	//  9          68719476736 0.4319874
  // 	// 10        1099511627776 0.2173495
  // 	// 11       17592186044416 0.110361
  // 	// 12      281474976710656 0.05478334
  // 	// 13     4503599627370496 0.02755045
  // 	// 14    72057594037927936 0.01351415
  // 	// 15  1152921504606846976 0.006807914
  // test cov rad
  {
    cout << "QuaternionMap Covrad" << endl;
    int NRES = 6;
    int NITER = 50000;
    NEST<4, Eigen::Quaterniond, QuaternionMap> nest;
    for (int r = 1; r <= NRES; ++r) {
      double maxdiff = 0, avgdiff = 0;
      for (int i = 0; i <= NITER; ++i) {
        Eigen::Quaterniond q(fabs(gauss(rng)), gauss(rng), gauss(rng),
                             gauss(rng));
        q.normalize();
        // cout << q << endl;
        Eigen::Quaterniond qcen = nest.set_and_get(nest.get_index(q, r), r);
        maxdiff = std::max(maxdiff, q.angularDistance(qcen));
        avgdiff += q.angularDistance(qcen);
      }
      avgdiff /= NITER;
      size_t count = 0;
      for (size_t i = 0; i < nest.size(r); ++i)
        if (nest.set_state(i, r)) ++count;
      double volfrac = (double)count * (maxdiff * maxdiff * maxdiff) * 4.0 /
                       3.0 * M_PI / 8.0 / M_PI / M_PI;
      double avgfrac = (double)count * (avgdiff * avgdiff * avgdiff) * 4.0 /
                       3.0 * M_PI / 8.0 / M_PI / M_PI;
      printf("%2i %16lu %10.5f %10.5f %10.5f %10.5f %10.5f\n", r, count,
             maxdiff * 180.0 / M_PI, avgdiff * 180.0 / M_PI, maxdiff / avgdiff,
             volfrac, avgfrac);
    }
  }

  // cout << "attempt to dump visualiztion" << endl;
  // {
  // 	typedef Eigen::Vector3d V;
  // 	NEST<4,Eigen::Quaterniond,QuaternionMap> nest;
  // 	// for(size_t r = 4; r < 5; ++r){
  // 	Eigen::Quaterniond qrand( gauss(rng), gauss(rng), gauss(rng), gauss(rng)
  // );
  // 	qrand.normalize();
  // 	Eigen::Vector3d X = qrand*Eigen::Vector3d(1,0,0);
  // 	Eigen::Vector3d Y = qrand*Eigen::Vector3d(0,1.1,0);
  // 	Eigen::Vector3d Z = qrand*Eigen::Vector3d(0,0,1.2);
  // 	size_t beg = 0;
  // 	while(!nest.set_state(beg,10)) beg =
  // std::max<size_t>(uniform(rng)*(nest.size(10)-1000),0);
  // 	for(size_t r = 1; r <= 9; ++r){
  // 		std::ofstream
  // out(("quatmaptest_"+boost::lexical_cast<std::string>(r)+".pdb").c_str());
  // 		size_t count = 0;
  // 		// cout << r << " " << nest.size(r) << " " << (beg>>(4*(10-r)))
  // <<
  // endl;
  // 		// continue;
  // 		for(size_t i = beg>>(4*(10-r)); i < nest.size(r); ++i){
  // 			if( nest.set_state(i,r) ){
  // 				++count;
  // 				if( count > 300) break;
  // 				Eigen::Quaterniond q = nest.value();
  // 				cout << r << " " << i << " " <<
  // q.coeffs().transpose()
  // <<
  // endl;
  // 				 V ximg = q.matrix() * X;
  // 				 V yimg = q.matrix() * Y;
  // 				 V zimg = q.matrix() * Z;;
  // 				 io::dump_pdb_atom(out,"O" ,61*ximg);
  // 				 io::dump_pdb_atom(out,"NI",61*yimg);
  // 				 io::dump_pdb_atom(out,"N" ,61*zimg);
  // 			}
  // 		}
  // 		out.close();
  // 	}
  // }
}
}
}
}
