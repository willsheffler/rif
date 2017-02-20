#include <gtest/gtest.h>

#include "numeric/euler_angles.hpp"

#include <Eigen/Geometry>

#include <boost/format.hpp>
#include <boost/lexical_cast.hpp>
#include <random>

#include "util/Timer.hpp"

namespace scheme {
namespace numeric {
namespace test {

using std::cout;
using std::endl;

TEST(euler_angles, test) {
  using namespace Eigen;
  std::mt19937 rng((unsigned int)time(0));
  std::normal_distribution<> gauss;
  std::uniform_real_distribution<> uniform;

  for (int i = 0; i < 10000; ++i) {
    Quaterniond q(fabs(gauss(rng)), gauss(rng), gauss(rng), gauss(rng));
    q.normalize();
    Matrix3d m = q.matrix();
    Vector3d euler;
    euler_angles(m, euler);
    // cout << euler.transpose() << endl;
    // cout << m << endl;
    // cout << endl;
    Matrix3d m2;
    from_euler_angles(euler, m2);
    double thresh = 0.0000001;
    if (euler[2] < 0.000001 || euler[2] > M_PI - 0.000001) thresh = 0.0002;
    if (euler[2] < 0.0000001 || euler[2] > M_PI - 0.0000001) thresh = 0.002;
    if (!m2.isApprox(m, thresh)) {
      cout << euler << endl;
      cout << m << endl;
      cout << m2 << endl;
    }
    ASSERT_TRUE(m2.isApprox(m, thresh));
  }
}

TEST(euler_angles, performance) {
  using namespace Eigen;
  std::mt19937 rng((unsigned int)time(0));
  std::normal_distribution<> gauss;
  std::uniform_real_distribution<> uniform;

  int NSAMP = 1 * 1000 * 1000;

  std::vector<Matrix3d> samp(NSAMP);
  std::vector<Vector3d> euler(NSAMP);

  for (int i = 0; i < NSAMP; ++i) {
    Quaterniond q(fabs(gauss(rng)), gauss(rng), gauss(rng), gauss(rng));
    q.normalize();
    Matrix3d m = q.matrix();
    samp[i] = m;
  }

  util::Timer<> t;
  for (int i = 0; i < NSAMP; ++i) {
    euler_angles(samp[i], euler[i]);
  }
  cout << "rate " << NSAMP / (double)t.elapsed() << endl;

  for (int i = 0; i < NSAMP; ++i) {
    // cout << euler.transpose() << endl;
    // cout << m << endl;
    // cout << endl;
    Matrix3d m2;
    from_euler_angles(euler[i], m2);
    double thresh = 0.0000001;
    if (euler[i][2] < 0.000001 || euler[i][2] > M_PI - 0.000001)
      thresh = 0.0002;
    if (euler[i][2] < 0.0000001 || euler[i][2] > M_PI - 0.0000001)
      thresh = 0.002;
    if (!m2.isApprox(samp[i], thresh)) {
      cout << euler[i] << endl;
      cout << samp[i] << endl;
      cout << m2 << endl;
    }
    ASSERT_TRUE(m2.isApprox(samp[i], thresh));
  }
}
}
}
}
