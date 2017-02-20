#include <gtest/gtest.h>

#include "bloom/bloom_filter.hpp"

#include <boost/foreach.hpp>
#include <random>
#include "util/Timer.hpp"

#include <sparsehash/dense_hash_map>

namespace scheme {
namespace util {
namespace bloom_test {

using std::cout;
using std::endl;

// avg num configurations: 1 1.55106
// avg num configurations: 2 2.40722
// avg num configurations: 3 3.73373
// avg num configurations: 4 5.792
// avg num configurations: 5 8.99201
// avg num configurations: 6 13.9436
// avg num configurations: 7 21.6277
// avg num configurations: 8 33.5358
// avg num configurations: 9 52.0504
// avg num configurations: 10 80.8304
// avg num configurations: 11 125.224
// avg num configurations: 12 194.104
// avg num configurations: 13 300.746
// avg num configurations: 14 467.801
// avg num configurations: 15 725.23
// avg num configurations: 16 1123.09
// avg num configurations: 17 1740.96
// avg num configurations: 18 2719.56
// avg num configurations: 19 4202.5

// TEST( TEMP_TEST, N_DOCK_DESIGN_GUESS ){
// 	double cdf[4];
//    	cdf[0] = 1.0;
//    	cdf[1] = ( 0.369952-0.0285942 ) / ( 1.0-0.0285942 );
//    	cdf[2] = ( 0.174289-0.0285942 ) / ( 1.0-0.0285942 );
//    	cdf[3] = ( 0.077201-0.0285942 ) / ( 1.0-0.0285942 );
//    	cout << cdf[3] << endl;
//    	cout << cdf[2] << endl;
//    	cout << cdf[1] << endl;
//    	cout << cdf[0] << endl;

// 	std::mt19937 rng((uint64_t)time(0));
// 	std::uniform_real_distribution<> uniform;

// 	int NSAMP = 10000000;

// 	for( int NRES = 1; NRES < 20; ++NRES ){
// 		double avg_tot_rot = 0;
// 		for( int i = 0; i < NSAMP; ++i ){
// 			std::vector<double> nrot(NRES);
// 			for( int j = 0; j < NRES; ++j ){
// 				double rand = uniform(rng);
// 				if     ( rand < cdf[3] ) nrot[j] = 4;
// 				else if( rand < cdf[2] ) nrot[j] = 3;
// 				else if( rand < cdf[1] ) nrot[j] = 2;
// 				else                     nrot[j] = 1;
// 			}
// 			int tot_rot = 1;
// 			for( int j = 0; j < nrot.size(); ++j ) tot_rot *=
// nrot[j];
// 			avg_tot_rot += tot_rot;
// 		}
// 		avg_tot_rot /= NSAMP;
// 		cout << "avg num configurations: " << NRES << " " << avg_tot_rot
// <<
// endl;
// 	}
// }

#ifdef SCHEME_BENCHMARK

TEST(bloom, bloom_filter_example) {
  uint64_t const NELEM = 100llu * 1000000llu;
  uint64_t const MAXIDX = NELEM * 1000000llu;
  uint64_t const NTEST = NELEM;

  bloom_parameters parameters;
  parameters.projected_element_count = NELEM;
  // parameters.optimal_parameters.table_size = 65536*8 * 512; // bits
  // parameters.optimal_parameters.number_of_hashes = 2;
  parameters.false_positive_probability = 0.5;
  parameters.random_seed = 0xA5A5A5A5;
  parameters.compute_optimal_parameters();

  std::cout << "Nelem: " << NELEM << " table_size: "
            << parameters.optimal_parameters.table_size / 8000000.0 << "M "
            << " Nhash: " << parameters.optimal_parameters.number_of_hashes
            << std::endl;

  std::mt19937 rng((uint64_t)0);
  std::uniform_int_distribution<int64_t> randindex(0, MAXIDX);

  ASSERT_FALSE(!parameters);

  bloom_filter filter(parameters);

  google::dense_hash_map<uint64_t, double> inserted;
  inserted.set_empty_key(std::numeric_limits<uint64_t>::max());
  for (uint64_t i = 0; i < NELEM; ++i) {
    uint64_t index = randindex(rng);
    filter.insert(index);
    inserted.insert(std::make_pair(index, (double)index));
  }

  // std::cout << "ninserted: " << inserted.size() << std::endl;
  // BOOST_FOREACH( uint64_t index, inserted ){
  // 	ASSERT_TRUE( filter.contains( index ) );
  // }

  {
    util::Timer<> timer;
    int nfalse = 0;
    for (uint64_t i = 0; i < NTEST; ++i) {
      nfalse += filter.contains(randindex(rng));
    }
    double t = timer.elapsed();
    cout << "false positive rate: " << (double)nfalse / NTEST << endl;
    cout << "BLOOM: " << NTEST / t / 1000000.0 << "M/s" << endl;
  }
  {
    util::Timer<> timer;
    int nfalse = 0;
    for (uint64_t i = 0; i < NTEST; ++i) {
      nfalse += inserted.find(randindex(rng)) != inserted.end();
    }
    double t = timer.elapsed();
    cout << "HASH: " << NTEST / t / 1000000.0 << "M/s" << endl;
  }
}

#endif
}
}
}
