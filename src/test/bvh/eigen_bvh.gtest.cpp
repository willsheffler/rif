#include <gtest/gtest.h>
#include <iostream>
#include <unsupported/Eigen/BVH>

#include "rif/eigen_types.hpp"
#include "rif/global_rng.hpp"

using namespace Eigen;
using namespace rif;

typedef AlignedBox<float, 3> Box3d;

namespace Eigen {
Box3d bounding_box(const V3f &v) { return Box3d(v, v); }
}

namespace rif_test_bvh_eigen_bvh {

struct PointPointMinimizer {
  typedef float Scalar;
  float minimumOnVolumeVolume(const Box3d &r1, const Box3d &r2) {
    ++calls;
    return r1.squaredExteriorDistance(r2);
  }
  float minimumOnVolumeObject(const Box3d &r, const V3f &v) {
    ++calls;
    return r.squaredExteriorDistance(v);
  }
  float minimumOnObjectVolume(const V3f &v, const Box3d &r) {
    ++calls;
    return r.squaredExteriorDistance(v);
  }
  float minimumOnObjectObject(const V3f &v1, const V3f &v2) {
    ++calls;
    return (v1 - v2).squaredNorm();
  }
  int calls = 0;
  X3f x = X3f::Identity();
};

TEST(eigen_bvh, example) {
  srand(global_rng()());
  typedef std::vector<V3f, aligned_allocator<V3f> > StdVectorOfVector3d;
  StdVectorOfVector3d redPoints, bluePoints;

  for (int i = 0; i < 1000; ++i) {
    redPoints.push_back(V3f::Random());
    bluePoints.push_back(V3f::Random() + V3f(0, 0, 0));
  }

  // brute force
  PointPointMinimizer minimizer;
  float brutemin = std::numeric_limits<float>::max();
  // brute force to find closest red-blue pair
  for (int i = 0; i < (int)redPoints.size(); ++i)
    for (int j = 0; j < (int)bluePoints.size(); ++j)
      brutemin = std::min(brutemin, minimizer.minimumOnObjectObject(
                                        redPoints[i], bluePoints[j]));
  int brutecalls = minimizer.calls;

  // bvh
  minimizer.calls = 0;
  KdBVH<float, 3, V3f> redTree(redPoints.begin(), redPoints.end()),
      blueTree(bluePoints.begin(), bluePoints.end());  // construct the trees
  float bvhmin = BVMinimize(redTree, blueTree, minimizer);
  int bvhcalls = minimizer.calls;

  std::cout << "             Brute force distance = " << sqrt(brutemin)
            << ", calls = " << brutecalls << std::endl;
  std::cout << "             BVH distance         = " << sqrt(bvhmin)
            << ", calls = " << bvhcalls << std::endl;
  std::cout << "             ratio: " << 1. * brutecalls / bvhcalls
            << std::endl;
}
}