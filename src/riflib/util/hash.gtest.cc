#ifdef CXX11

#include <gtest/gtest.h>

#include <random>
#include "riflib/util/Timer.hh"
#include <boost/foreach.hpp>

#include <sparsehash/dense_hash_map>
#include <sparsehash/sparse_hash_map>

#include "riflib/util/SimpleArray.hh"

#include <unordered_map>
#include <boost/unordered_map.hpp>
#include <map>

namespace scheme { namespace util { namespace test_hash {

using std::cout;
using std::endl;

template<class K, class V, class GetMap, size_t SEG=0>
struct SegmentedMap {
	typedef typename GetMap::template apply< K, util::SimpleArray<1<<SEG,V,true> >::type MAP;
	typedef typename MAP::const_iterator const_iterator;

	MAP map;

	SegmentedMap() { map.set_empty_key(std::numeric_limits<K>::max()); }

	const_iterator find(K const & k) const {
		return map.find(k>>SEG);
	}

	const_iterator end() const {
		return map.end();
	}

	void insert(std::pair<K,V> const & v){
		map.insert(std::make_pair(v.first>>SEG,v.second));
	}

	V const & operator[](K const & k) const {
		return map.find(k>>SEG)->second.operator[]( k % (1<<SEG) );
	}

	V & operator[](K const & k) {
		return map[k>>SEG][ k % (1<<SEG) ];
	}
};


template<class Map>
void test_map( Map * hp, int64_t MAXIDX, int64_t NSAMP ){
	Map & h(*hp);

	int64_t NROW=1;
	NSAMP /= NROW;

	// std::mt19937 rng((unsigned int64_t)time(0));
	std::mt19937 rng((uint64_t)0);	
	std::uniform_int_distribution<int64_t> randindex(0,MAXIDX);
	// h.resize(NFILL/2);

	util::Timer<> t;
	size_t count = 0;
	for(size_t i = 0; i < NSAMP; ++i){
		size_t ri = randindex(rng);
		for(size_t j = 0; j < NROW; ++j){
			typename Map::const_iterator iter = h.find(ri+j);
			count += iter==h.end() ? 0.0 : iter->second[0];
		}
	}

	if( count == 0 ) cout << "ERROR!!" << endl;
	cout << "rate "  << NSAMP << " " << NROW << " " << MAXIDX << " " << t.elapsed_nano()/NSAMP/NROW << "ns, nonsense: " << (double)count/NSAMP << endl;
}


template<class Map>
void fill_map(Map & h, int64_t MAXIDX, int64_t sparsity=100ll ){
	int64_t NFILL = MAXIDX/sparsity;

	std::mt19937 rng((uint64_t)0);	
	std::uniform_int_distribution<int64_t> randindex(0,MAXIDX);
	// h.resize(NFILL/2);

	for(int64_t i = 0; i < NFILL; ++i) h[randindex(rng)] = i;
	cout << "done fill " 
	     << MAXIDX << " " 
	     << (double)h.size()/1000000 << "M  "
	     << (double)h.bucket_count()/1000000000 * sizeof(typename Map::value_type) << "GB  " 
	     << (double)h.size() / h.bucket_count()
	     << endl;
}


template<class Map>
void fill_and_test_map(Map * hp){
	Map & h(*hp);
	int64_t MAXIDX = 3000ll*1000ll*1000ll;
	int64_t NFILL = MAXIDX/100ll;
	int64_t NROW = 1;
	int64_t NSAMP = 10ll*1000ll*1000ll / NROW;

	std::mt19937 rng((uint64_t)0);	
	std::uniform_int_distribution<int64_t> randindex(0,MAXIDX);
	// h.resize(NFILL/2);

	for(int64_t i = 0; i < NFILL; ++i) h[randindex(rng)] = i;
	cout << "done fill " 
	     << MAXIDX << " "
	     << (double)h.size()/1000000 << "M  "
	     << (double)h.bucket_count()/1000000000 * sizeof(typename Map::value_type) << "GB  " 
	     << (double)h.size() / h.bucket_count()
	     << endl;	     	     

	std::mt19937 rng2((uint64_t)0);	

	util::Timer<> t;
	size_t count = 0;
	for(size_t i = 0; i < NSAMP; ++i){
		size_t ri = randindex(rng2);
		for(size_t j = 0; j < NROW; ++j){
			typename Map::const_iterator iter = h.find(ri+j);
			count += iter==h.end() ? 0.0 : iter->second[0];
		}
	}
	cout << "rate " << NSAMP << " " << NROW << " " << MAXIDX << " " << t.elapsed_nano()/NSAMP/NROW << "ns, nonsense: " << (double)count / NSAMP << endl;
}

// 100M 1/10 nrow 1
// dense rate 1.04204e+07 -6.16073
// std   rate 4.29674e+06 -18.0023
// boost rate 5.60315e+06 4.98979

// 100M 1/10 nrow 10
// dense rate 2.15978e+07 176.247
// boost rate 6.85248e+06 200.992

struct GoogleDense { template<class K, class V> struct apply { typedef google::dense_hash_map<K,V> type; }; };

// TEST(test_hash, sparse_vs_dense){
// 	// SegmentedMap<uint64_t,util::SimpleArray<8,double> ,GoogleDense,2> m;

// 	google::dense_hash_map<uint64_t,util::SimpleArray<8,double> > d;
// 	d.set_empty_key(std::numeric_limits<uint64_t>::max());

// 	// google::sparse_hash_map<uint64_t,util::SimpleArray<8,double> > s;

// 	// fill_and_test_map(n,"google_dense",m,"segment_gdh ");
// 	cout << "====================== DENSE ======================" << endl;

// 	int64_t MAXIDX = 3000ll*1000ll*1000ll;
// 	int64_t NSAMP = 10ll*1000ll*1000ll;
// 	fill_map( d, MAXIDX, 100 );
// 	test_map( &d, MAXIDX, NSAMP );
// 	d.clear();

// 	fill_and_test_map(&d);
// 	d.clear();

// 	// cout << "====================== SPARSE =====================" << endl;
// 	// fill_and_test_map(s); s.clear();

// }



template<class MAP1, class MAP2>
void test_2map(
	MAP1 & m, std::string lm,
	MAP2 & n, std::string ln
){
	int64_t MAXIDX = 10*1000*1000;
	int64_t NFILL = MAXIDX/10;
	int64_t NROW = 1;
	int64_t NSAMP = 100*1000*1000 / NROW;

	// std::mt19937 rng((unsigned int64_t)time(0));
	std::mt19937 rng((uint64_t)0);	
	std::uniform_int_distribution<int64_t> randindex(0,MAXIDX);

	for(int64_t i = 0; i < NFILL; ++i){
		size_t ri = randindex(rng);
		m[ri] = 1;
		n[ri] = 1;
	}

	std::vector<size_t> ridx;
	for(size_t i = 0; i < NSAMP; ++i)
		ridx.push_back( randindex(rng) );

	util::Timer<> tm;
	size_t mcount = 0;
	BOOST_FOREACH(size_t ri,ridx){
		for(size_t j = 0; j < NROW; ++j){
			if( m.find(ri+j) != m.end() ) mcount += m[ri+j];
		}
	}
	cout << lm << " rate " << tm.elapsed_nano()/NSAMP*NROW << "ns" << endl;

	util::Timer<> tn;
	size_t ncount = 0;
	BOOST_FOREACH(size_t ri,ridx){
		for(size_t j = 0; j < NROW; ++j){
			if( n.find(ri+j) != n.end() ) ncount += n[ri+j];
		}
	}
	cout << ln << " rate " << tn.elapsed_nano()/NSAMP*NROW << "ns" << endl;

	cout << mcount << endl;
	ASSERT_EQ( mcount, ncount );
}



// TEST(test_hash,std_unordered_map){
// 	std::unordered_map<int64_t,int64_t> h;
// 	fill_and_test_map(h);
// }

// TEST(test_hash,boost_unordered_map){
// 	boost::unordered_map<int64_t,int64_t> h;
// 	fill_and_test_map(h);
// }

// TEST(test_hash,std_map){
// 	std::map<int64_t,int64_t> h;
// 	fill_and_test_map(h);
// }


}}}

#endif

