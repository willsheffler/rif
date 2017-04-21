#include <gtest/gtest.h>
#include <iostream>
#include <random>

#include <rif/geom/bvh.hpp>
#include "rif/geom/primitive.hpp"
#include "rif/global_rng.hpp"
#include "rif/numeric/rand_xform.hpp"
#include "util/Timer.hpp"

using namespace Eigen;
using namespace rif;
using namespace rif::geom;
using namespace rif::geom::bvh;

using F = double;
using Timer = util::Timer;

namespace Eigen {
auto bounding_vol(V3d v) { return Sphered(v); }
}

namespace rif_test_bvh_eigen_bvh {

struct PPMin {
  typedef F Scalar;
  using Xform = X3<F>;
  PPMin(Xform x = Xform::Identity()) : bXa(x) {}
  Scalar minimumOnVolumeVolume(Sphere<Scalar> r1, Sphere<Scalar> r2) {
    ++calls;
    return r1.signdis(bXa * r2);
  }
  Scalar minimumOnVolumeObject(Sphere<Scalar> r, V3d v) {
    ++calls;
    return r.signdis(bXa * v);
  }
  Scalar minimumOnObjectVolume(V3d v, Sphere<Scalar> r) {
    ++calls;
    return (bXa * r).signdis(v);
  }
  Scalar minimumOnObjectObject(V3d v1, V3d v2) {
    ++calls;
    return (v1 - bXa * v2).norm();
  }
  int calls = 0;
  Xform bXa = Xform::Identity();
  void reset() { calls = 0; }
};

TEST(eigen_bvh, test_min) {
  typedef std::vector<V3d, aligned_allocator<V3d> > StdVectorOfVector3d;
  StdVectorOfVector3d ptsA, ptsB;
  std::uniform_real_distribution<> r(0, 1);
  std::mt19937& g(global_rng());
  for (double dx = 0.91; dx < 1.1; dx += 0.02) {
    StdVectorOfVector3d ptsA, ptsB;
    for (int i = 0; i < 100; ++i) {
      ptsA.push_back(V3d(r(g), r(g), r(g)));
      ptsB.push_back(V3d(r(g), r(g), r(g)) + V3d(dx, 0, 0));
    }

    // brute force
    PPMin minimizer;
    auto tbrute = Timer("tb");
    F brutemin = std::numeric_limits<F>::max();
    // brute force to find closest red-blue pair
    for (int i = 0; i < (int)ptsA.size(); ++i)
      for (int j = 0; j < (int)ptsB.size(); ++j)
        brutemin = std::min(brutemin,
                            minimizer.minimumOnObjectObject(ptsA[i], ptsB[j]));
    tbrute.stop();
    int brutecalls = minimizer.calls;

    // bvh
    // move Pa by random X, set bXa in minimizer
    auto X = rif::numeric::rand_xform(F(999));
    for (auto& p : ptsA) p = X * p;
    minimizer.bXa = X;

    minimizer.reset();
    auto tcreate = Timer("tc");
    WelzlBVH<double, V3d> bvhA(ptsA.begin(), ptsA.end()),
        bvhB(ptsB.begin(), ptsB.end());  // construct the trees
    tcreate.stop();
    auto tvbh = Timer("tbvh");
    F bvhmin = BVMinimize(bvhA, bvhB, minimizer);
    tvbh.stop();
    int bvhcalls = minimizer.calls;

    ASSERT_FLOAT_EQ(brutemin, bvhmin);

    float ratio = 1. * brutecalls / bvhcalls;
    std::cout << "    Brute/BVH " << dx << " " << ratio << " " << brutemin
              << " " << bvhmin << " " << brutecalls << " " << bvhcalls << " "
              << tbrute << " " << tcreate << std::endl;
  }
}

struct PPIsect {
  using Scalar = F;
  using Xform = X3<F>;
  PPIsect(F r, Xform x = Xform::Identity()) : radius(r), bXa(x) {}
  bool intersectVolumeVolume(Sphere<Scalar> r1, Sphere<Scalar> r2) {
    ++calls;
    return r1.signdis(bXa * r2) < radius;
  }
  bool intersectVolumeObject(Sphere<Scalar> r, V3d v) {
    ++calls;
    return r.signdis(bXa * v) < radius;
  }
  bool intersectObjectVolume(V3d v, Sphere<Scalar> r) {
    ++calls;
    return (bXa * r).signdis(v) < radius;
  }
  bool intersectObjectObject(V3d v1, V3d v2) {
    ++calls;
    bool isect = (v1 - bXa * v2).norm() < radius;
    result |= isect;
    return isect;
  }
  void reset() {
    calls = 0;
    result = false;
  }
  int calls = 0;
  F radius = 0.0;
  bool result = false;
  Xform bXa = Xform::Identity();
};

TEST(eigen_bvh, test_isect) {
  typedef std::vector<V3d, aligned_allocator<V3d> > StdVectorOfVector3d;
  std::uniform_real_distribution<> r(0, 1);
  std::mt19937& g(global_rng());
  double avg_ratio = 0.0;
  int niter = 0;
  for (double dx = 0.001 + 0.95; dx < 1.05; dx += 0.005) {
    ++niter;

    StdVectorOfVector3d ptsA, ptsB;
    for (int i = 0; i < 100; ++i) {
      ptsA.push_back(V3d(r(g), r(g), r(g)));
      ptsB.push_back(V3d(r(g), r(g), r(g)) + V3d(dx, 0, 0));
    }
    PPIsect query(0.1);

    // brute force
    bool bruteisect = false;
    // brute force to find closest red-blue pair
    auto tbrute = Timer("tb");
    for (int i = 0; i < (int)ptsA.size(); ++i) {
      for (int j = 0; j < (int)ptsB.size(); ++j) {
        if (query.intersectObjectObject(ptsA[i], ptsB[j])) {
          bruteisect = true;
          break;
        }
      }
      if (bruteisect) break;
    }
    int brutecalls = query.calls;
    tbrute.stop();

    // bvh

    query.reset();

    // WTF??? this does not cause the test to fail...
    auto X = rif::numeric::rand_xform(F(999));
    for (auto& p : ptsA) p = X * p;
    query.bXa = X;  // commenting this out should fail

    auto tcreate = Timer("tc");
    WelzlBVH<double, V3d> bvhA(ptsA.begin(), ptsA.end()),
        bvhB(ptsB.begin(), ptsB.end());
    tcreate.stop();
    // std::cout << bvhA.vols[0] << std::endl;
    // std::cout << bvhB.vols[0] << std::endl;

    auto tbvh = Timer("tbvh");
    BVIntersect(bvhA, bvhB, query);
    tbvh.stop();
    bool bvhisect = query.result;
    int bvhcalls = query.calls;

    ASSERT_EQ(bruteisect, bvhisect);

    float ratio = 1. * brutecalls / bvhcalls;
    avg_ratio += ratio;
    std::cout << "    Brute/BVH " << dx << " " << ratio << " " << bruteisect
              << " " << bvhisect << " " << brutecalls << " " << bvhcalls << " "
              << tbrute << " " << tcreate << std::endl;
  }
  avg_ratio /= niter;
  std::cout << "avg Brute/BVH " << avg_ratio << std::endl;
}
}
