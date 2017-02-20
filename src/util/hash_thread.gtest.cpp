#ifdef CXX11
#ifdef SCHEME_BENCHMARK

#include <gtest/gtest.h>

#include <boost/foreach.hpp>
#include <boost/multi_array.hpp>
#include <random>
#include "util/Timer.hpp"

#include <sparsehash/dense_hash_map>
#include <sparsehash/sparse_hash_map>

#include "util/SimpleArray.hpp"

#include <boost/unordered_map.hpp>
#include <map>
#include <unordered_map>

#include <mutex>
#include <thread>

namespace scheme {
namespace util {
namespace test_hash_thread {

using std::cout;
using std::endl;

template <class K, class V, class GetMap, size_t SEG = 0>
struct SegmentedMap {
  typedef typename GetMap::template apply<
      K, util::SimpleArray<1 << SEG, V, true> >::type MAP;
  typedef typename MAP::const_iterator const_iterator;

  MAP map;

  SegmentedMap() { map.set_empty_key(std::numeric_limits<K>::max()); }

  const_iterator find(K const& k) const { return map.find(k >> SEG); }

  const_iterator end() const { return map.end(); }

  void insert(std::pair<K, V> const& v) {
    map.insert(std::make_pair(v.first >> SEG, v.second));
  }

  V const& operator[](K const& k) const {
    return map.find(k >> SEG)->second.operator[](k % (1 << SEG));
  }

  V& operator[](K const& k) { return map[k >> SEG][k % (1 << SEG)]; }
};

template <class Map>
void fill_map(Map& h, int64_t MAXIDX, int64_t sparsity = 100ll) {
  int64_t NFILL = MAXIDX / sparsity;

  std::mt19937 rng((uint64_t)0);
  std::uniform_int_distribution<int64_t> randindex(0, MAXIDX);
  // h.resize(NFILL/2);

  util::Timer<> t;
  for (int64_t i = 0; i < NFILL; ++i) h[randindex(rng)] = i;
  cout << "done fill " << (double)h.size() / 1000000 << "M  entries, "
       << (double)h.bucket_count() / 1000000000 *
              sizeof(typename Map::value_type)
       << "GB  total size, " << (double)h.size() / h.bucket_count() << " load, "
       << t.elapsed_nano() / NFILL << "ns" << endl;
}

template <class Map>
void test_map(Map* hp, double* runtime, int64_t MAXIDX, int64_t NITER,
              std::mt19937* rng0) {
  Map const& h(*hp);

  int64_t NROW = 1;
  NITER /= NROW;

  static std::mutex m;

  std::uniform_int_distribution<int64_t> randindex(0, MAXIDX);
  std::mt19937 rng(randindex(*rng0));

  // h.resize(NFILL/2);

  util::Timer<> t;
  size_t count = 0;
  for (size_t i = 0; i < NITER; ++i) {
    size_t ri = randindex(rng);
    for (size_t j = 0; j < NROW; ++j) {
      typename Map::const_iterator iter = h.find(ri + j);
      count += iter == h.end() ? 0.0 : iter->second[0];
    }
  }
  double const time = t.elapsed_nano();

  m.lock();
  *runtime += time;
  if (count == 0) cout << "ERROR!!" << endl;
  // cout << "rate " << t.elapsed_nano()/NITER/NROW << "ns, nonsense: " <<
  // (double)count << endl;
  m.unlock();
}

struct GoogleDense {
  template <class K, class V>
  struct apply {
    typedef google::dense_hash_map<K, V> type;
  };
};

#ifdef SCHEME_BENCHMARK

TEST(test_hash_thread, google_hash_thread) {
  int64_t MAXIDX = 2000ll * 1000ll * 1000ll;
  int64_t NSAMP_TOT = 100ll * 1000ll * 1000ll;
  int maxNthread = 16;

  std::mt19937 rng((unsigned int)time(0));
  typedef google::dense_hash_map<uint64_t, util::SimpleArray<8, double> > D;
  D d;
  d.set_empty_key(std::numeric_limits<uint64_t>::max());

  cout << "====================== HASH_TEST ======================" << endl;
  fill_map(d, MAXIDX, 100);

  for (int Nthread = 1; Nthread <= maxNthread; ++Nthread) {
    int64_t NSAMP = NSAMP_TOT / Nthread;
    // test_map(d); d.clear();
    double runtime = 0.0;
    std::vector<std::thread> t;
    for (int i = 0; i < Nthread; ++i)
      t.push_back(std::thread(test_map<D>, &d, &runtime, MAXIDX, NSAMP, &rng));
    for (int i = 0; i < Nthread; ++i) t[i].join();
    runtime /= Nthread;
    printf(
        "nthread: %5i  %7.3f ns / lookup  %7.3f s runtime  %7.3f M lookup/sec  "
        "%7.3f M lookup/sec/thread \n",
        Nthread, runtime / NSAMP_TOT, runtime / 1000000000.0,
        NSAMP_TOT / runtime * 1000.0, NSAMP_TOT / runtime / Nthread * 1000.0);
    std::cout.flush();
  }
  d.clear();
}

#endif
// on lappy 4x haswell 2.8
// nthread:     1  118.039 ns / lookup   11.804 s runtime    8.472 M lookup/sec
// 8.472 M lookup/sec/thread
// nthread:     2   56.901 ns / lookup    5.690 s runtime   17.575 M lookup/sec
// 8.787 M lookup/sec/thread
// nthread:     3   36.720 ns / lookup    3.672 s runtime   27.233 M lookup/sec
// 9.078 M lookup/sec/thread
// nthread:     4   32.093 ns / lookup    3.209 s runtime   31.160 M lookup/sec
// 7.790 M lookup/sec/thread
// nthread:     5   29.967 ns / lookup    2.997 s runtime   33.370 M lookup/sec
// 6.674 M lookup/sec/thread
// nthread:     6   28.283 ns / lookup    2.828 s runtime   35.357 M lookup/sec
// 5.893 M lookup/sec/thread
// nthread:     7   26.248 ns / lookup    2.625 s runtime   38.098 M lookup/sec
// 5.443 M lookup/sec/thread
// nthread:     8   25.216 ns / lookup    2.522 s runtime   39.657 M lookup/sec
// 4.957 M lookup/sec/thread
// nthread:     9   24.903 ns / lookup    2.490 s runtime   40.156 M lookup/sec
// 4.462 M lookup/sec/thread
// nthread:    10   24.859 ns / lookup    2.486 s runtime   40.227 M lookup/sec
// 4.023 M lookup/sec/thread
// nthread:    11   24.995 ns / lookup    2.500 s runtime   40.008 M lookup/sec
// 3.637 M lookup/sec/thread
// nthread:    12   25.353 ns / lookup    2.535 s runtime   39.443 M lookup/sec
// 3.287 M lookup/sec/thread
// nthread:    13   25.427 ns / lookup    2.543 s runtime   39.328 M lookup/sec
// 3.025 M lookup/sec/thread
// nthread:    14   25.538 ns / lookup    2.554 s runtime   39.158 M lookup/sec
// 2.797 M lookup/sec/thread
// nthread:    15   25.138 ns / lookup    2.514 s runtime   39.781 M lookup/sec
// 2.652 M lookup/sec/thread
// nthread:    16   24.468 ns / lookup    2.447 s runtime   40.870 M lookup/sec
// 2.554 M lookup/sec/thread

void test_array(float const* const h, double* runtime, size_t N, int64_t NITER,
                std::mt19937* rng) {
  std::uniform_int_distribution<int64_t> randindex(0, N);
  // h.resize(NFILL/2);

  util::Timer<> t;
  float count = 0;
  for (size_t i = 0; i < NITER; ++i) {
    size_t ri = randindex(*rng);
    count += *(h + ri);
  }
  double const time = t.elapsed_nano();

  static std::mutex m;
  m.lock();
  *runtime += time;
  if (count == 0) cout << "ERROR!!" << endl;
  // cout << "rate " << t.elapsed_nano()/NITER/NROW << "ns, nonsense: " <<
  // (double)count << endl;
  m.unlock();
}

// real	8m41.960s
// 6m22.756s

TEST(test_hash_thread, simple_array_thread) {
  int maxNthread = 8;
  int64_t SIZE = 1000 * 1000 * 1000;
  // int64_t SIZE = 100*100*10*5;
  int64_t NSAMP_TOT = 100ll * 1000ll * 1000ll;

  std::vector<float> data(SIZE);
  for (int64_t i = 0; i < SIZE; ++i) {
    data[i] = (float)i;
  }
  cout << "done fill " << (double)SIZE * 4.0 / 1000000000.0 << "GB" << endl;

  std::mt19937 rng((unsigned int)time(0));

  // test_map(n,"google_dense",m,"segment_gdh ");
  cout << "====================== ARRAY_TEST ======================" << endl;

  // double rt=0;
  // test_map<D>( &d, &rt, MAXIDX, NSAMP_TOT );
  // cout << "main thread: " << rt << "ns / lookup" << endl;

  for (int Nthread = 1; Nthread <= maxNthread; ++Nthread) {
    int64_t NSAMP = NSAMP_TOT / Nthread;
    // test_map(d); d.clear();
    double runtime = 0.0;
    std::vector<std::thread> t;
    for (int i = 0; i < Nthread; ++i)
      t.push_back(
          std::thread(test_array, &data[0], &runtime, SIZE, NSAMP, &rng));
    for (int i = 0; i < Nthread; ++i) t[i].join();
    runtime /= Nthread;
    printf(
        "nthread: %5i, %7.3fns / lookup, %7.3fs runtime, %7.3fM lookup/sec, "
        "%7.3fM lookup/sec/thread \n",
        Nthread, runtime / NSAMP_TOT, runtime / 1000000000.0,
        NSAMP_TOT / runtime * 1000.0, NSAMP_TOT / runtime / Nthread * 1000.0);
    std::cout.flush();
  }
  data.clear();
}

void test_multiarray(boost::multi_array<float, 1> const* hp, double* runtime,
                     int64_t NITER, std::mt19937* rng) {
  std::uniform_int_distribution<int64_t> randindex(0, hp->shape()[0]);
  // h.resize(NFILL/2);

  util::Timer<> t;
  float count = 0;
  for (size_t i = 0; i < NITER; ++i) {
    size_t ri = randindex(*rng);
    count += (*hp)[ri];
  }
  double const time = t.elapsed_nano();

  static std::mutex m;
  m.lock();
  *runtime += time;
  if (count == 0) cout << "ERROR!!" << endl;
  // cout << "rate " << t.elapsed_nano()/NITER/NROW << "ns, nonsense: " <<
  // (double)count << endl;
  m.unlock();
}

// real	8m41.960s
// 6m22.756s

#ifdef SCHEME_BENCHMARK

TEST(test_hash_thread, multi_array_thread) {
  int maxNthread = 32;
  int64_t SIZE = 1000 * 1000 * 1000;
  // int64_t SIZE = 100*100*10*5;
  int64_t NSAMP_TOT = 100ll * 1000ll * 1000ll;

  boost::multi_array<float, 1> data(boost::extents[SIZE]);

  for (int64_t i = 0; i < SIZE; ++i) {
    data[i] = (float)i;
  }
  cout << "done fill " << (double)SIZE * 4.0 / 1000000000.0 << "GB" << endl;

  std::mt19937 rng((unsigned int)time(0));

  // test_map(n,"google_dense",m,"segment_gdh ");
  cout << "====================== MULTI_ARRAY_TEST ======================"
       << endl;

  // double rt=0;
  // test_map<D>( &d, &rt, MAXIDX, NSAMP_TOT );
  // cout << "main thread: " << rt << "ns / lookup" << endl;

  for (int Nthread = 1; Nthread <= maxNthread; ++Nthread) {
    int64_t NSAMP = NSAMP_TOT / Nthread;
    // test_map(d); d.clear();
    double runtime = 0.0;
    std::vector<std::thread> t;
    for (int i = 0; i < Nthread; ++i)
      t.push_back(std::thread(test_multiarray, &data, &runtime, NSAMP, &rng));
    for (int i = 0; i < Nthread; ++i) t[i].join();
    runtime /= Nthread;
    printf(
        "nthread: %5i, %7.3fns / lookup, %7.3fs runtime, %7.3fM lookup/sec, "
        "%7.3fM lookup/sec/thread \n",
        Nthread, runtime / NSAMP_TOT, runtime / 1000000000.0,
        NSAMP_TOT / runtime * 1000.0, NSAMP_TOT / runtime / Nthread * 1000.0);
    std::cout.flush();
  }
}

#endif

// nthread:     1,  84.195ns / lookup,   8.420s runtime,  11.877M lookup/sec,
// 11.877M lookup/sec/thread
// nthread:     2,  67.559ns / lookup,   6.756s runtime,  14.802M lookup/sec,
// 7.401M lookup/sec/thread
// nthread:     3,  51.433ns / lookup,   5.143s runtime,  19.443M lookup/sec,
// 6.481M lookup/sec/thread
// nthread:     4,  43.524ns / lookup,   4.352s runtime,  22.976M lookup/sec,
// 5.744M lookup/sec/thread
// nthread:     5,  37.039ns / lookup,   3.704s runtime,  26.998M lookup/sec,
// 5.400M lookup/sec/thread
// nthread:     6,  33.654ns / lookup,   3.365s runtime,  29.714M lookup/sec,
// 4.952M lookup/sec/thread
// nthread:     7,  31.637ns / lookup,   3.164s runtime,  31.609M lookup/sec,
// 4.516M lookup/sec/thread
// nthread:     8,  31.232ns / lookup,   3.123s runtime,  32.018M lookup/sec,
// 4.002M lookup/sec/thread

// 	// cout << "====================== SPARSE =====================" <<
// endl;
// 	// test_map(s); s.clear();
// }

// ====================== dense_hash_map ======================
// done fill 39.8005M  9.66368GB  0.296537583.638ns
// main thread: 117.639
//                 dense       std    sparse
//                 9.7gb     3.8gb     3.0gb
// runtime     1 116.507   327.728   275.698
// runtime     2  47.328   158.355   115.406
// runtime     3  33.946   109.843    77.206
// runtime     4  24.162    86.198    58.460
// runtime     5  18.175    68.707    45.108
// runtime     6  16.432    57.417    39.264
// runtime     7  16.162    49.047    34.603
// runtime     8  23.386    46.593    34.271
// runtime     9  20.376    39.326    29.788
// runtime    10  15.526    32.694    27.142
// runtime    11  14.557    29.697    26.838
// runtime    12  14.538    29.900    26.953
// runtime    13  15.157    30.562    26.904
// runtime    14  14.894    32.081    25.822
// runtime    15  15.061    30.378    26.667
// runtime    16  15.734    31.775    27.109
// runtime    17  14.636    30.269    26.191
// runtime    18  14.716    29.658    25.827
// runtime    19  15.201    31.496    25.571
// runtime    20  14.980    34.030    25.837
// runtime    21  14.604    29.368    25.897
// runtime    22  14.157    28.862    25.810
// runtime    23  15.036    30.119    25.423
// runtime    24  14.653    31.193    26.734
// runtime    25  14.665    32.219    29.968
// runtime    26  14.626    31.588    30.335
// runtime    27  14.553    32.113    25.477
// runtime    28  14.699    34.224    26.479
// runtime    29  14.010    29.427    30.412
// runtime    30  14.786    29.614    27.794
// runtime    31  14.423    30.900    28.131
// runtime    32  14.913    32.057    27.559

// ====================== std::unordered_map ======================
// done fill 39.8005M  3.79296GB  0.755515664.705ns
// main thread: 315.627
// runtime     1 327.728
// runtime     2 158.355
// runtime     3 109.843
// runtime     4  86.198
// runtime     5  68.707
// runtime     6  57.417
// runtime     7  49.047
// runtime     8  46.593
// runtime     9  39.326
// runtime    10  32.694
// runtime    11  29.697
// runtime    12  29.900
// runtime    13  30.562
// runtime    14  32.081
// runtime    15  30.378
// runtime    16  31.775
// runtime    17  30.269
// runtime    18  29.658
// runtime    19  31.496
// runtime    20  34.030
// runtime    21  29.368
// runtime    22  28.862
// runtime    23  30.119
// runtime    24  31.193
// runtime    25  32.219
// runtime    26  31.588
// runtime    27  32.113
// runtime    28  34.224
// runtime    29  29.427
// runtime    30  29.614
// runtime    31  30.900
// runtime    32  32.057

// ====== sparsehash ACTUAL MEM ~3.0gb
// done fill 39.8005M  4.83184GB  0.5930741179.81ns
// main thread: 279.367
// runtime     1 275.698
// runtime     2 115.406
// runtime     3  77.206
// runtime     4  58.460
// runtime     5  45.108
// runtime     6  39.264
// runtime     7  34.603
// runtime     8  34.271
// runtime     9  29.788
// runtime    10  27.142
// runtime    11  26.838
// runtime    12  26.953
// runtime    13  26.904
// runtime    14  25.822
// runtime    15  26.667
// runtime    16  27.109
// runtime    17  26.191
// runtime    18  25.827
// runtime    19  25.571
// runtime    20  25.837
// runtime    21  25.897
// runtime    22  25.810
// runtime    23  25.423
// runtime    24  26.734
// runtime    25  29.968
// runtime    26  30.335
// runtime    27  25.477
// runtime    28  26.479
// runtime    29  30.412
// runtime    30  27.794
// runtime    31  28.131
// runtime    32  27.559

// dense_hash_map
// runtime     1 186.508 298.519
// runtime     2  85.484 141.353
// runtime     3  56.238  97.312
// runtime     4  43.064  74.743
// runtime     5  36.117  59.752
// runtime     6  33.066  53.263
// runtime     7  30.003  47.826
// runtime     8  32.299  45.008

// std::unordered_map
// runtime     1 298.519
// runtime     2 141.353
// runtime     3  97.312
// runtime     4  74.743
// runtime     5  59.752
// runtime     6  53.263
// runtime     7  47.826
// runtime     8  45.008
// TEST(test_hash,std_unordered_map){
// 	std::unordered_map<int64_t,int64_t> h;
// 	test_map(h);
// }

// TEST(test_hash,boost_unordered_map){
// 	boost::unordered_map<int64_t,int64_t> h;
// 	test_map(h);
// }

// TEST(test_hash,std_map){
// 	std::map<int64_t,int64_t> h;
// 	test_map(h);
// }
}
}
}

#endif
#endif
