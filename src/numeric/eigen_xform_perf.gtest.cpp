// #define EIGEN_DONT_VECTORIZE 1

#include <gtest/gtest.h>

#include "numeric/rand_xform.hpp"

#include <Eigen/Geometry>

#include <random>
#include "util/Timer.hpp"

namespace rif {
namespace numeric {
namespace pref_test {

using std::cout;
using std::endl;

typedef Eigen::Transform<double, 3, Eigen::AffineCompact> Xform;
// typedef Eigen::Affine3d Xform;

// TEST( XformMap, basic_test ){
template <class Xform>
void test_xform_perf() {
  int NSAMP = 10 * 1000;
#ifdef SCHEME_BENCHMARK
  NSAMP = 1 * 1000 * 1000;
#endif
  std::mt19937 rng;
  Xform x, sum = Xform::Identity();
  rand_xform(rng, x);

  double mintime = 9e9, tottime = 0;
  for (int iter = 0; iter < 10; ++iter) {
    util::Timer<> t;
    for (int i = 0; i < NSAMP; ++i) {
      sum = sum * x;
    }
    double time = t.elapsed_nano();
    mintime = std::min(time, mintime);
    tottime += time;
  }
  printf("runtime %7.3fns nonsense: %7.3f tot: %f\n", mintime / NSAMP,
         sum.translation()[0], tottime);
}

TEST(xform_perf, preformance) {
#ifdef EIGEN_VECTORIZE
  cout << "EIGEN_VECTORIZE is set" << endl;
#else
  cout << "EIGEN_VECTORIZE is NOT set" << endl;
#endif
  cout << "AffineCompact d ";
  test_xform_perf<Eigen::Transform<double, 3, Eigen::AffineCompact>>();
  cout << "Affine        d ";
  test_xform_perf<Eigen::Transform<double, 3, Eigen::Affine>>();
  cout << "AffineCompact f ";
  test_xform_perf<Eigen::Transform<float, 3, Eigen::AffineCompact>>();
  cout << "Affine        f ";
  test_xform_perf<Eigen::Transform<float, 3, Eigen::Affine>>();
  // cout << "XformHash_bt24_Cubic_Zorder"; test_xform_perf<
  // XformHash_bt24_Cubic_Zorder >();
}
}
}
}
