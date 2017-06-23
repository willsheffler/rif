#include <gtest/gtest.h>

#include "io/dump_pdb_atom.hpp"
#include "nest/pmap/TetracontoctachoronMap.hpp"
#include "numeric/lattice.hpp"

#include <boost/foreach.hpp>
#include <random>

#include <sparsehash/dense_hash_set>

#include <Eigen/Geometry>

namespace rif {
namespace numeric {
namespace test {

using std::cout;
using std::endl;

// QuaternionMap Covrad
//  1                8  119.36447   61.22694    1.94954    3.83748    0.51790
//  2              128   59.04495   27.31148    2.16191    7.43172    0.73549
//  3             1040   27.85890   13.60802    2.04724    6.34244    0.73918
//  4             8480   13.59714    6.82994    1.99081    6.01270    0.76204
//  5            68416    6.91475    3.40731    2.02939    6.37995    0.76335
//  6           548496    3.44178    1.70610    2.01734    6.30746    0.76828
//  7          4392160    1.73778    0.85268    2.03801    6.50123    0.76802

void test_orientatin_coverage_4d(size_t Nside, int NSAMP) {
  typedef double F;
  typedef uint64_t S;
  typedef util::SimpleArray<4, F> V;
  typedef util::SimpleArray<4, S> I;
  std::mt19937 rng((unsigned int)time(0));
  std::normal_distribution<> rnorm;
  std::uniform_real_distribution<> runif;

  // Cubic<4,F,S> bcc(I(Nside),V(-1),V(1.0));
  // BCC<4,F,S> bcc(I(Nside+2),V(-1.0-2.0/Nside),V(1.0+2.0/Nside));
  BCC<4, F, S> bcc(I(Nside + 1), V(-1.0 - 2.0 / Nside), V(1.0));

  google::dense_hash_set<size_t> idx_seen;
  idx_seen.set_empty_key(999999999999);
  F maxdiff = 0, avgdiff = 0;
  for (int i = 0; i < NSAMP; ++i) {
    V sample(rnorm(rng), rnorm(rng), rnorm(rng), rnorm(rng));
    sample /= sample.norm();
    if (sample[0] < 0) sample = -sample;
    ASSERT_DOUBLE_EQ(sample.norm(), 1.0);
    size_t index = bcc[sample];
    idx_seen.insert(index);
    V center = bcc[index];
    Eigen::Quaterniond samp(sample[0], sample[1], sample[2], sample[3]);
    Eigen::Quaterniond cen(center[0], center[1], center[2], center[3]);
    cen.normalize();
    maxdiff = std::max(maxdiff, samp.angularDistance(cen));
    avgdiff += samp.angularDistance(cen);
  }
  avgdiff /= NSAMP;
  size_t count = idx_seen.size();
  double xcov = (double)NSAMP / count;
  if (xcov < 40) {
    cout << "volfrac will be low if xcov < 50 or so... (not enough samples)"
         << endl;
    cout << "idxfrac: " << (double)count / (Nside * Nside * Nside * Nside)
         << " xcov " << xcov << endl;
  }
  double volfrac = (double)count * (maxdiff * maxdiff * maxdiff) * 4.0 / 3.0 *
                   M_PI / 8.0 / M_PI / M_PI;
  double avgfrac = (double)count * (avgdiff * avgdiff * avgdiff) * 4.0 / 3.0 *
                   M_PI / 8.0 / M_PI / M_PI;
  // cout << "RESOL         COUNT   CovRad      AvgRad     mx/avg   VolFrac
  // AvgFrac" << endl;
  printf("%2lu %16lu %10.5f %10.5f %10.5f %10.5f %10.5f\n", Nside, count,
         maxdiff * 180.0 / M_PI, avgdiff * 180.0 / M_PI, maxdiff / avgdiff,
         volfrac, avgfrac);

  // ASSERT_LE( volfrac , 3.5 );
}

void test_orientatin_coverage_3d_bt24(size_t Nside, int NSAMP) {
  typedef double F;
  typedef uint64_t S;
  typedef util::SimpleArray<4, F> V;
  typedef util::SimpleArray<4, S> I;
  typedef util::SimpleArray<3, F> V3;
  typedef util::SimpleArray<3, S> I3;
  std::mt19937 rng((unsigned int)time(0));
  // std::mt19937 rng((unsigned int)0);
  std::normal_distribution<> rnorm;

  // BCC<3,F,S> bcc(I3(Nside),V3( 0.0 ),V3( 1.0));
  // BCC<3,F,S> bcc(I3(3*Nside),V3(-1.0 ),V3( 2.0));

  BCC<3, F, S> bcc(I3(Nside + 1), V3(-1.0 / Nside), V3(1.0));
  // Cubic<3,F,S> bcc(I3(Nside),V3(0.0),V3(1.0));

  nest::pmap::TetracontoctachoronMap<> map;
  typedef nest::pmap::TetracontoctachoronMap<3, V>::Params Params;
  // cout << bcc.lower_ << endl;
  // cout << bcc.width_ << endl;
  // cout << bcc.lower_+bcc.width_*Nside << endl;

  google::dense_hash_set<size_t> idx_seen;
  idx_seen.set_empty_key(999999999999);
  F maxdiff = 0, avgdiff = 0;
  for (int i = 0; i < NSAMP; ++i) {
    V sample(rnorm(rng), rnorm(rng), rnorm(rng), rnorm(rng));
    sample /= sample.norm();
    if (sample[0] < 0) sample = -sample;
    ASSERT_DOUBLE_EQ(sample.norm(), 1.0);

    Eigen::Quaterniond q;
    for (int i = 0; i < 4; ++i) q.coeffs()[i] = sample[i];

    Params params;
    uint64_t cell_index;
    map.value_to_params(q.matrix(), 0, params, cell_index);
    size_t index_bcc = bcc[params];
    assert(index_bcc < bcc.size());
    idx_seen.insert(index_bcc);
    // size_t index = index_bcc + bcc.size()*cell_index;
    Params pcenter = bcc[index_bcc];
    // if( !( pcenter[0] >= 0.0 && pcenter[0] <= 1.0 ) ) cout << "bad cell cen "
    // << Nside << " " << pcenter << endl;
    // if( !( pcenter[1] >= 0.0 && pcenter[1] <= 1.0 ) ) cout << "bad cell cen "
    // << Nside << " " << pcenter << endl;
    // if( !( pcenter[2] >= 0.0 && pcenter[2] <= 1.0 ) ) cout << "bad cell cen "
    // << Nside << " " << pcenter << endl;
    // assert( pcenter[0] >= 0.0 && pcenter[0] <= 1.0 );
    // assert( pcenter[1] >= 0.0 && pcenter[1] <= 1.0 );
    // assert( pcenter[2] >= 0.0 && pcenter[2] <= 1.0 );
    Eigen::Matrix3d mcenter;
    // cout << pcenter << " " << cell_index << endl;
    bool params_to_value_returns_success =
        map.params_to_value(pcenter, cell_index, 0, mcenter);
    assert(params_to_value_returns_success);

    Eigen::Quaterniond qcenter(mcenter);
    V center;
    for (int i = 0; i < 4; ++i) center[i] = qcenter.coeffs()[i];

    Eigen::Quaterniond samp(sample[0], sample[1], sample[2], sample[3]);
    Eigen::Quaterniond cen(center[0], center[1], center[2], center[3]);
    cen.normalize();
    maxdiff = std::max(maxdiff, samp.angularDistance(cen));
    avgdiff += samp.angularDistance(cen);
  }
  avgdiff /= NSAMP;
  // size_t count = bcc.size() * 24;
  size_t count = 24 * idx_seen.size();
  double xcov = (double)NSAMP / count * 24;
  if (xcov < 50) {
    cout << "volfrac will be low if xcov < 50 or so... (not enough samples)"
         << endl;
    cout << "idxfrac: " << (double)count / (Nside * Nside * Nside) << " xcov "
         << xcov << endl;
  }
  double volfrac = (double)count * (maxdiff * maxdiff * maxdiff) * 4.0 / 3.0 *
                   M_PI / 8.0 / M_PI / M_PI;
  double avgfrac = (double)count * (avgdiff * avgdiff * avgdiff) * 4.0 / 3.0 *
                   M_PI / 8.0 / M_PI / M_PI;
  // cout << "RESOL         COUNT   CovRad      AvgRad     mx/avg   VolFrac
  // AvgFrac" << endl;
  printf("%2lu %16lu %10.5f %10.5f %10.5f %10.5f %10.5f \n", Nside, count,
         maxdiff * 180.0 / M_PI, avgdiff * 180.0 / M_PI, maxdiff / avgdiff,
         volfrac, avgfrac);

  // ASSERT_LE( volfrac , 3.5 );
}

#ifdef SCHEME_BENCHMARK
int NSAMP = 600 * 1000;
int NSIDE = 12;
#else
int NSAMP = 30 * 1000;
int NSIDE = 6;
#endif

TEST(bcc_lattice, orientatin_coverage_4d) {
  cout << "RESOL         COUNT     CovRad     AvgRad     mx/avg    VolFrac    "
          "AvgFrac"
       << endl;
  for (int i = 3; i < 3 * NSIDE - 3; i += 3) {
    test_orientatin_coverage_4d(i, 6 * NSAMP);
    std::cout.flush();
  }
  //  test_orientatin_coverage_4d(   8 , NSAMP );
  //  test_orientatin_coverage_4d(  16 , NSAMP );
  //  test_orientatin_coverage_4d(  32 , NSAMP );
  //  // test_orientatin_coverage_4d(  64 , NSAMP );
  //  // test_orientatin_coverage_4d( 128 , NSAMP );
}

TEST(bcc_lattice, orientatin_coverage_3d_bt24) {
  std::mt19937 rng((unsigned int)time(0));
  std::uniform_real_distribution<> runif;
  cout << "RESOL         COUNT     CovRad     AvgRad     mx/avg    VolFrac    "
          "AvgFrac"
       << endl;
  for (int i = 1; i < NSIDE; ++i) {
    test_orientatin_coverage_3d_bt24(i, NSAMP);
    std::cout.flush();
  }
  // test_orientatin_coverage_3d_bt24(  4 , NSAMP );
  // test_orientatin_coverage_3d_bt24(  8 , NSAMP );
  // test_orientatin_coverage_3d_bt24( 16 , NSAMP );
  // test_orientatin_coverage_3d_bt24( 32 , NSAMP );
  // test_orientatin_coverage_3d_bt24( 64 , NSAMP );
}

// CUBIC:
// RESOL         COUNT     CovRad     AvgRad     mx/avg    VolFrac    AvgFrac
//  1               24   62.74152   40.73375    1.54028    1.67189    0.45752
//  2              192   39.31569   20.63676    1.90513    3.29102    0.47595
//  3              648   26.76186   13.82667    1.93552    3.50312    0.48312
//  4             1536   20.19943   10.38852    1.94440    3.57059    0.48572
//  5             3000   16.18746    8.31584    1.94658    3.58912    0.48660
//  6             5184   13.55390    6.92925    1.95604    3.64074    0.48647
//  7             8232   11.62521    5.94226    1.95636    3.64787    0.48718
//  8            12288   10.10018    5.19960    1.94249    3.57108    0.48722
//  9            17496    9.04248    4.62274    1.95609    3.64865    0.48749
// 10            23880    8.11903    4.15979    1.95179    3.60478    0.48482

// BCC
// RESOL         COUNT     CovRad     AvgRad     mx/avg    VolFrac    AvgFrac
//  1              216   49.63837   31.72439    1.56468    7.45139    1.94520
//  2              840   25.97057   16.08611    1.61447    4.15006    0.98620
//  3             1992   17.41470   10.70403    1.62693    2.96735    0.68907
//  4             4344   13.12570    8.04887    1.63075    2.77069    0.63889
//  5             7992   10.50987    6.44961    1.62954    2.61685    0.60477
//  6            12648    8.74000    5.37186    1.62700    2.38170    0.55300
//  7            19752    7.50751    4.60429    1.63055    2.35738    0.54379
//  8            28824    6.55751    4.03079    1.62686    2.29246    0.53242
//  9            40536    5.84932    3.58299    1.63252    2.28817    0.52591
// 10            53832    5.27287    3.22479    1.63510    2.22594    0.50919

// CUBIC  1               24   62.74152   40.73375    1.54028    1.67189 0.45752
// BCC   1              216   49.63837   31.72439    1.56468    7.45139 1.94520
// CUBIC  2              192   39.31569   20.63676    1.90513    3.29102 0.47595
// CUBIC  3              648   26.76186   13.82667    1.93552    3.50312 0.48312
// BCC   2              840   25.97057   16.08611    1.61447    4.15006 0.98620
// CUBIC  4             1536   20.19943   10.38852    1.94440    3.57059 0.48572
// BCC   3             1992   17.41470   10.70403    1.62693    2.96735 0.68907
// CUBIC  5             3000   16.18746    8.31584    1.94658    3.58912 0.48660
// CUBIC  6             5184   13.55390    6.92925    1.95604    3.64074 0.48647
// BCC   4             4344   13.12570    8.04887    1.63075    2.77069 0.63889
// CUBIC  7             8232   11.62521    5.94226    1.95636    3.64787 0.48718
// BCC   5             7992   10.50987    6.44961    1.62954    2.61685 0.60477
// CUBIC  8            12288   10.10018    5.19960    1.94249    3.57108 0.48722
// CUBIC  9            17496    9.04248    4.62274    1.95609    3.64865 0.48749
// BCC   6            12648    8.74000    5.37186    1.62700    2.38170 0.55300
// CUBIC 10            23880    8.11903    4.15979    1.95179    3.60478 0.48482
// BCC   7            19752    7.50751    4.60429    1.63055    2.35738 0.54379
// BCC   8            28824    6.55751    4.03079    1.62686    2.29246 0.53242
// BCC   9            40536    5.84932    3.58299    1.63252    2.28817 0.52591
// BCC  10            53832    5.27287    3.22479    1.63510    2.22594 0.50919

// conclusion, use 4D for covering > ~15°, and bt24 3D for < ~15°
// or figure out how to better optimize 48-cell based lookup, esp. at low grid
// number
// RESOL         COUNT   CovRad      AvgRad     mx/avg   VolFrac      AvgFrac
//  2                1        nan        nan        nan        nan        nan
//  3               33   74.28860   48.00784    1.54743    3.81603    1.02987
//  4               85   54.29327   30.41631    1.78500    3.83697    0.67464
//  5              169   43.45715   25.11014    1.73066    3.91202    0.75468
//  6              325   33.33101   19.38110    1.71977    3.39437    0.66734
//  7              576   29.51199   16.11477    1.83136    4.17589    0.67987
//  8              877   24.50840   14.06207    1.74287    3.64146    0.68783
//  9             1310   21.81570   12.10485    1.80223    3.83627    0.65536
// 10             1915   19.12114   10.54955    1.81251    3.77608    0.63416
// 11             2749   16.93394    9.54341    1.77441    3.76513    0.67393
// 12             3283   15.63030    8.66045    1.80479    3.53593    0.60148
// 13             4451   14.02718    7.85187    1.78647    3.46497    0.60773
// 14             5581   13.09922    7.26046    1.80419    3.53817    0.60247
// 15             7241   11.91341    6.71108    1.77519    3.45333    0.61731
// 16             8755   11.14681    6.20759    1.79567    3.42010    0.59068
// 17            10921   10.36544    5.81305    1.78313    3.43049    0.60507
// 18            13046    9.84143    5.49671    1.79042    3.50738    0.61111
// 19            15131    9.12415    5.15155    1.77114    3.24172    0.58346
}
}
}
