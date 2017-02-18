#include <gtest/gtest.h>

#include "types.hpp"
#include "nest/MultiNest.hpp"

#include <boost/any.hpp>
#include <iterator>     // std::ostream_iterator

namespace scheme { namespace nest { namespace test {

using std::cout;
using std::endl;


TEST( MultiNest, boost_any ){
	using boost::any;
	using boost::any_cast;

	double d = 123;
	double *dp = &d;
	shared_ptr<double> dsp = make_shared<double>();

	any ad = d;
	ASSERT_EQ( d, any_cast<double>(ad) );
	ASSERT_EQ( d, *any_cast<double>(&ad) );
	ASSERT_EQ( nullptr, any_cast<float>(&ad) );

	any adp = dp;

	d = 1.2345; ASSERT_EQ( 1.2345, *any_cast<double*>(adp) );
	d = 3.2345; ASSERT_EQ( 3.2345, *any_cast<double*>(adp) );
	d = 3.2345; ASSERT_EQ( 3.2345, **any_cast<double*>(&adp) );
	d = 3.2345; ASSERT_EQ( nullptr, any_cast<float*>(&adp) );

	any adsp = dsp;
	*dsp = 2.34; ASSERT_EQ( 2.34, *any_cast<shared_ptr<double> >(adsp) );
	*dsp = 3.34; ASSERT_EQ( 3.34, **any_cast<shared_ptr<double> >(&adsp) );
	*dsp = 2.34; ASSERT_EQ( nullptr, any_cast<double*>(&adsp) );
 }

TEST( MultiNest, expand_index ){
	std::ostream_iterator<uint64_t> out_it (std::cout,", ");
	MultiNest<uint64_t,uint64_t>::Indices indices;

	MultiNest<uint64_t,uint64_t> multi;
	ASSERT_EQ( multi.virtual_dim(), 0 );

	indices.clear();
	multi.expand_index( 0, indices );
	ASSERT_EQ( indices.size(), 0 );

	multi.add_nest( make_shared<NEST<1> >() );
	ASSERT_EQ( multi.virtual_dim(), 1 );

	for(int i = 0; i < 10; ++i){
		indices.clear();
		multi.expand_index( i, indices );
		ASSERT_EQ( indices.size(), 1 );
		ASSERT_EQ( indices[0], i );
	}

	multi.add_nest( make_shared<NEST<2> >() );
	ASSERT_EQ( multi.virtual_dim(), 3 );

	std::cout << endl;
	int sum1=0, sum2=0, sum3=0;
	for(int i = 0; i < 64; ++i){
		indices.clear();
		multi.expand_index( i, indices );
		ASSERT_EQ( indices.size(), 3 );
		sum1 += indices[0];
		sum2 += indices[1];
		sum3 += indices[2];
		// ASSERT_EQ( indices[0], i );
		// std::cout << "RESULT: ";
		// std::copy ( indices.begin(), indices.end(), out_it );
		// std::cout << endl;
	}
	ASSERT_EQ( sum1, 64*3/2 );
	ASSERT_EQ( sum1, 64*3/2 );
	ASSERT_EQ( sum1, 64*3/2 );
 }

TEST( MultiNest, get_state_single_cell ){
	std::ostream_iterator<uint64_t> out_it (std::cout,", ");
	MultiNest<uint64_t,uint64_t>::Indices indices;

	MultiNest<uint64_t,uint64_t> multi;
	multi.add_nest( make_shared<NEST<1> >() );
	multi.add_nest( make_shared<NEST<2> >() );
	multi.add_nest( make_shared<NEST<3> >() );
	ASSERT_EQ(multi.size(0),1);

	typedef util::SimpleArray<1,double> A1; A1 val1;
	typedef util::SimpleArray<2,double> A2; A2 val2;
	typedef util::SimpleArray<3,double> A3; A3 val3;
	std::vector<boost::any> anys(3);
	anys[0] = &val1;
	anys[1] = &val2;
	anys[2] = &val3;

	multi.get_states( 0, 0, anys ); ASSERT_EQ( val1, A1( 0.5 ) ); ASSERT_EQ( val2, A2( 0.5, 0.5 ) ); ASSERT_EQ( val3, A3( 0.5, 0.5, 0.5 ) );

	multi.get_states( 0, 1, anys ); ASSERT_EQ( val1, A1( 0.25) ); ASSERT_EQ( val2, A2( 0.25, 0.25 ) ); ASSERT_EQ( val3, A3( 0.25, 0.25, 0.25 ) );
	multi.get_states( 1, 1, anys ); ASSERT_EQ( val1, A1( 0.75) ); ASSERT_EQ( val2, A2( 0.25, 0.25 ) ); ASSERT_EQ( val3, A3( 0.25, 0.25, 0.25 ) );
	multi.get_states( 2, 1, anys ); ASSERT_EQ( val1, A1( 0.25) ); ASSERT_EQ( val2, A2( 0.75, 0.25 ) ); ASSERT_EQ( val3, A3( 0.25, 0.25, 0.25 ) );
	multi.get_states( 5, 1, anys ); ASSERT_EQ( val1, A1( 0.75) ); ASSERT_EQ( val2, A2( 0.25, 0.75 ) ); ASSERT_EQ( val3, A3( 0.25, 0.25, 0.25 ) );

	multi.get_states( 0, 2, anys ); ASSERT_EQ( val1, A1( 0.125 ) ); ASSERT_EQ( val2, A2( 0.125, 0.125 ) ); ASSERT_EQ( val3, A3( 0.125, 0.125, 0.125 ) );
	multi.get_states( 1, 2, anys ); ASSERT_EQ( val1, A1( 0.375 ) ); ASSERT_EQ( val2, A2( 0.125, 0.125 ) ); ASSERT_EQ( val3, A3( 0.125, 0.125, 0.125 ) );
	multi.get_states( 2, 2, anys ); ASSERT_EQ( val1, A1( 0.125 ) ); ASSERT_EQ( val2, A2( 0.375, 0.125 ) ); ASSERT_EQ( val3, A3( 0.125, 0.125, 0.125 ) );
	multi.get_states( 5, 2, anys ); ASSERT_EQ( val1, A1( 0.375 ) ); ASSERT_EQ( val2, A2( 0.125, 0.375 ) ); ASSERT_EQ( val3, A3( 0.125, 0.125, 0.125 ) );

	multi.get_states( (1<<6)+0, 2, anys ); ASSERT_EQ( val1, A1( 0.625 ) ); ASSERT_EQ( val2, A2( 0.125, 0.125 ) ); ASSERT_EQ( val3, A3( 0.125, 0.125, 0.125 ) );
	multi.get_states( (1<<6)+1, 2, anys ); ASSERT_EQ( val1, A1( 0.875 ) ); ASSERT_EQ( val2, A2( 0.125, 0.125 ) ); ASSERT_EQ( val3, A3( 0.125, 0.125, 0.125 ) );
	multi.get_states( (1<<6)+2, 2, anys ); ASSERT_EQ( val1, A1( 0.625 ) ); ASSERT_EQ( val2, A2( 0.375, 0.125 ) ); ASSERT_EQ( val3, A3( 0.125, 0.125, 0.125 ) );
	multi.get_states( (1<<6)+5, 2, anys ); ASSERT_EQ( val1, A1( 0.875 ) ); ASSERT_EQ( val2, A2( 0.125, 0.375 ) ); ASSERT_EQ( val3, A3( 0.125, 0.125, 0.125 ) );
 }

TEST( MultiNest, get_index_single_cell ){
	std::ostream_iterator<uint64_t> out_it (std::cout,", ");
	MultiNest<uint64_t,uint64_t>::Indices indices;

	MultiNest<uint64_t,uint64_t> multi;
	multi.add_nest( make_shared<NEST<1> >() );
	multi.add_nest( make_shared<NEST<2> >() );
	multi.add_nest( make_shared<NEST<3> >() );
	ASSERT_EQ(multi.size(0),1);

	typedef util::SimpleArray<1,double> A1; A1 val1;
	typedef util::SimpleArray<2,double> A2; A2 val2;
	typedef util::SimpleArray<3,double> A3; A3 val3;
	std::vector<boost::any> anys(3);
	anys[0] = &val1;
	anys[1] = &val2;
	anys[2] = &val3;

	for(int resl = 0; resl < 4; ++resl){
		for(int index = 0; index < multi.size(resl); ++index){
			multi.get_states( index, resl, anys );
			// cout << "VALS: (" << val1 << ") (" << val2 << ") (" << val3 << ")"<<endl;
			uint64_t test_index = multi.virtual_get_index( &anys, resl );
			ASSERT_EQ( index, test_index );
		}
	}

 }

TEST( MultiNest, get_state_two_cells ){
	std::ostream_iterator<uint64_t> out_it (std::cout,", ");
	MultiNest<uint64_t,uint64_t>::Indices indices;

	MultiNest<uint64_t,uint64_t> multi;
	multi.add_nest( make_shared<NEST<1> >(2) );
	multi.add_nest( make_shared<NEST<2> >() );
	multi.add_nest( make_shared<NEST<3> >() );
	ASSERT_EQ(multi.size(0),2);

	typedef util::SimpleArray<1,double> A1; A1 val1;
	typedef util::SimpleArray<2,double> A2; A2 val2;
	typedef util::SimpleArray<3,double> A3; A3 val3;
	std::vector<boost::any> anys(3);
	anys[0] = &val1;
	anys[1] = &val2;
	anys[2] = &val3;

	int resl = 0;
	int h = (1<<multi.dim()*resl);

	multi.get_states( 0, resl, anys ); ASSERT_EQ( val1, A1( 0.5 ) ); ASSERT_EQ( val2, A2( 0.5, 0.5  ) ); ASSERT_EQ( val3, A3( 0.5, 0.5, 0.5 ) );

	multi.get_states( 0+h, resl, anys ); ASSERT_EQ( val1, A1( 1.5 ) ); ASSERT_EQ( val2, A2( 0.5, 0.5  ) ); ASSERT_EQ( val3, A3( 0.5, 0.5, 0.5 ) );
	resl = 1;
	h = (1<<multi.dim()*resl);

	multi.get_states( 0, 1, anys ); ASSERT_EQ( val1, A1( 0.25 ) ); ASSERT_EQ( val2, A2( 0.25, 0.25 ) ); ASSERT_EQ( val3, A3( 0.25, 0.25, 0.25 ) );
	multi.get_states( 1, 1, anys ); ASSERT_EQ( val1, A1( 0.75 ) ); ASSERT_EQ( val2, A2( 0.25, 0.25 ) ); ASSERT_EQ( val3, A3( 0.25, 0.25, 0.25 ) );
	multi.get_states( 2, 1, anys ); ASSERT_EQ( val1, A1( 0.25 ) ); ASSERT_EQ( val2, A2( 0.75, 0.25 ) ); ASSERT_EQ( val3, A3( 0.25, 0.25, 0.25 ) );
	multi.get_states( 5, 1, anys ); ASSERT_EQ( val1, A1( 0.75 ) ); ASSERT_EQ( val2, A2( 0.25, 0.75 ) ); ASSERT_EQ( val3, A3( 0.25, 0.25, 0.25 ) );

	multi.get_states( h+0, 1, anys ); ASSERT_EQ( val1, A1( 1.25 ) ); ASSERT_EQ( val2, A2( 0.25, 0.25 ) ); ASSERT_EQ( val3, A3( 0.25, 0.25, 0.25 ) );
	multi.get_states( h+1, 1, anys ); ASSERT_EQ( val1, A1( 1.75 ) ); ASSERT_EQ( val2, A2( 0.25, 0.25 ) ); ASSERT_EQ( val3, A3( 0.25, 0.25, 0.25 ) );
	multi.get_states( h+2, 1, anys ); ASSERT_EQ( val1, A1( 1.25 ) ); ASSERT_EQ( val2, A2( 0.75, 0.25 ) ); ASSERT_EQ( val3, A3( 0.25, 0.25, 0.25 ) );
	multi.get_states( h+5, 1, anys ); ASSERT_EQ( val1, A1( 1.75 ) ); ASSERT_EQ( val2, A2( 0.25, 0.75 ) ); ASSERT_EQ( val3, A3( 0.25, 0.25, 0.25 ) );

	multi.get_states( 0, 2, anys ); ASSERT_EQ( val1, A1( 0.125  ) ); ASSERT_EQ( val2, A2( 0.125, 0.125 ) ); ASSERT_EQ( val3, A3( 0.125, 0.125, 0.125 ) );
	multi.get_states( 1, 2, anys ); ASSERT_EQ( val1, A1( 0.375  ) ); ASSERT_EQ( val2, A2( 0.125, 0.125 ) ); ASSERT_EQ( val3, A3( 0.125, 0.125, 0.125 ) );
	multi.get_states( 2, 2, anys ); ASSERT_EQ( val1, A1( 0.125  ) ); ASSERT_EQ( val2, A2( 0.375, 0.125 ) ); ASSERT_EQ( val3, A3( 0.125, 0.125, 0.125 ) );
	multi.get_states( 5, 2, anys ); ASSERT_EQ( val1, A1( 0.375  ) ); ASSERT_EQ( val2, A2( 0.125, 0.375 ) ); ASSERT_EQ( val3, A3( 0.125, 0.125, 0.125 ) );

	multi.get_states( (1<<6)+0, 2, anys ); ASSERT_EQ( val1, A1( 0.625  ) ); ASSERT_EQ( val2, A2( 0.125, 0.125 ) ); ASSERT_EQ( val3, A3( 0.125, 0.125, 0.125 ) );
	multi.get_states( (1<<6)+1, 2, anys ); ASSERT_EQ( val1, A1( 0.875  ) ); ASSERT_EQ( val2, A2( 0.125, 0.125 ) ); ASSERT_EQ( val3, A3( 0.125, 0.125, 0.125 ) );
	multi.get_states( (1<<6)+2, 2, anys ); ASSERT_EQ( val1, A1( 0.625  ) ); ASSERT_EQ( val2, A2( 0.375, 0.125 ) ); ASSERT_EQ( val3, A3( 0.125, 0.125, 0.125 ) );
	multi.get_states( (1<<6)+5, 2, anys ); ASSERT_EQ( val1, A1( 0.875  ) ); ASSERT_EQ( val2, A2( 0.125, 0.375 ) ); ASSERT_EQ( val3, A3( 0.125, 0.125, 0.125 ) );
 }

TEST( MultiNest, get_index_two_cells ){
	std::ostream_iterator<uint64_t> out_it (std::cout,", ");
	MultiNest<uint64_t,uint64_t>::Indices indices;

	MultiNest<uint64_t,uint64_t> multi;
	multi.add_nest( make_shared<NEST<1> >(2) );
	multi.add_nest( make_shared<NEST<2> >() );
	multi.add_nest( make_shared<NEST<3> >() );
	ASSERT_EQ(multi.size(0),2);

	typedef util::SimpleArray<1,double> A1; A1 val1;
	typedef util::SimpleArray<2,double> A2; A2 val2;
	typedef util::SimpleArray<3,double> A3; A3 val3;
	std::vector<boost::any> anys(3);
	anys[0] = &val1;
	anys[1] = &val2;
	anys[2] = &val3;

	for(int resl = 0; resl < 3; ++resl){
		for(int index = 0; index < multi.size(resl); ++index){
			multi.get_states( index, resl, anys );
			// cout << "VALS: (" << val1 << ") (" << val2 << ") (" << val3 << ")"<<endl;
			uint64_t test_index = multi.virtual_get_index( &anys, resl );
			ASSERT_EQ( index, test_index );
		}
	}

 }

TEST( MultiNest, get_state_ncell_handling ){
	std::ostream_iterator<uint64_t> out_it (std::cout,", ");
	MultiNest<uint64_t,uint64_t>::Indices indices;

	MultiNest<uint64_t,uint64_t> multi;
	multi.add_nest( make_shared<NEST<1> >(4) );
	multi.add_nest( make_shared<NEST<2> >(3) );
	multi.add_nest( make_shared<NEST<3> >(2) );
	ASSERT_EQ(multi.size(0),24);

	typedef util::SimpleArray<1,double> A1; A1 val1;
	typedef util::SimpleArray<2,double> A2; A2 val2;
	typedef util::SimpleArray<3,double> A3; A3 val3;
	std::vector<boost::any> anys(3);
	anys[0] = &val1;
	anys[1] = &val2;
	anys[2] = &val3;

	multi.get_states(  0, 0, anys ); ASSERT_EQ( val1, A1( 0.5 ) ); ASSERT_EQ( val2, A2( 0.5, 0.5 ) ); ASSERT_EQ( val3, A3( 0.5, 0.5, 0.5 ) );
	multi.get_states(  1, 0, anys ); ASSERT_EQ( val1, A1( 1.5 ) ); ASSERT_EQ( val2, A2( 0.5, 0.5 ) ); ASSERT_EQ( val3, A3( 0.5, 0.5, 0.5 ) );
	multi.get_states(  2, 0, anys ); ASSERT_EQ( val1, A1( 2.5 ) ); ASSERT_EQ( val2, A2( 0.5, 0.5 ) ); ASSERT_EQ( val3, A3( 0.5, 0.5, 0.5 ) );
	multi.get_states(  3, 0, anys ); ASSERT_EQ( val1, A1( 3.5 ) ); ASSERT_EQ( val2, A2( 0.5, 0.5 ) ); ASSERT_EQ( val3, A3( 0.5, 0.5, 0.5 ) );
	multi.get_states(  4, 0, anys ); ASSERT_EQ( val1, A1( 0.5 ) ); ASSERT_EQ( val2, A2( 1.5, 0.5 ) ); ASSERT_EQ( val3, A3( 0.5, 0.5, 0.5 ) );
	multi.get_states(  5, 0, anys ); ASSERT_EQ( val1, A1( 1.5 ) ); ASSERT_EQ( val2, A2( 1.5, 0.5 ) ); ASSERT_EQ( val3, A3( 0.5, 0.5, 0.5 ) );
	multi.get_states(  6, 0, anys ); ASSERT_EQ( val1, A1( 2.5 ) ); ASSERT_EQ( val2, A2( 1.5, 0.5 ) ); ASSERT_EQ( val3, A3( 0.5, 0.5, 0.5 ) );
	multi.get_states(  7, 0, anys ); ASSERT_EQ( val1, A1( 3.5 ) ); ASSERT_EQ( val2, A2( 1.5, 0.5 ) ); ASSERT_EQ( val3, A3( 0.5, 0.5, 0.5 ) );
	multi.get_states(  8, 0, anys ); ASSERT_EQ( val1, A1( 0.5 ) ); ASSERT_EQ( val2, A2( 2.5, 0.5 ) ); ASSERT_EQ( val3, A3( 0.5, 0.5, 0.5 ) );
	multi.get_states(  9, 0, anys ); ASSERT_EQ( val1, A1( 1.5 ) ); ASSERT_EQ( val2, A2( 2.5, 0.5 ) ); ASSERT_EQ( val3, A3( 0.5, 0.5, 0.5 ) );
	multi.get_states( 10, 0, anys ); ASSERT_EQ( val1, A1( 2.5 ) ); ASSERT_EQ( val2, A2( 2.5, 0.5 ) ); ASSERT_EQ( val3, A3( 0.5, 0.5, 0.5 ) );
	multi.get_states( 11, 0, anys ); ASSERT_EQ( val1, A1( 3.5 ) ); ASSERT_EQ( val2, A2( 2.5, 0.5 ) ); ASSERT_EQ( val3, A3( 0.5, 0.5, 0.5 ) );
	multi.get_states( 12, 0, anys ); ASSERT_EQ( val1, A1( 0.5 ) ); ASSERT_EQ( val2, A2( 0.5, 0.5 ) ); ASSERT_EQ( val3, A3( 1.5, 0.5, 0.5 ) );
	multi.get_states( 13, 0, anys ); ASSERT_EQ( val1, A1( 1.5 ) ); ASSERT_EQ( val2, A2( 0.5, 0.5 ) ); ASSERT_EQ( val3, A3( 1.5, 0.5, 0.5 ) );
	multi.get_states( 14, 0, anys ); ASSERT_EQ( val1, A1( 2.5 ) ); ASSERT_EQ( val2, A2( 0.5, 0.5 ) ); ASSERT_EQ( val3, A3( 1.5, 0.5, 0.5 ) );
	multi.get_states( 15, 0, anys ); ASSERT_EQ( val1, A1( 3.5 ) ); ASSERT_EQ( val2, A2( 0.5, 0.5 ) ); ASSERT_EQ( val3, A3( 1.5, 0.5, 0.5 ) );
	multi.get_states( 16, 0, anys ); ASSERT_EQ( val1, A1( 0.5 ) ); ASSERT_EQ( val2, A2( 1.5, 0.5 ) ); ASSERT_EQ( val3, A3( 1.5, 0.5, 0.5 ) );
	multi.get_states( 17, 0, anys ); ASSERT_EQ( val1, A1( 1.5 ) ); ASSERT_EQ( val2, A2( 1.5, 0.5 ) ); ASSERT_EQ( val3, A3( 1.5, 0.5, 0.5 ) );
	multi.get_states( 18, 0, anys ); ASSERT_EQ( val1, A1( 2.5 ) ); ASSERT_EQ( val2, A2( 1.5, 0.5 ) ); ASSERT_EQ( val3, A3( 1.5, 0.5, 0.5 ) );
	multi.get_states( 19, 0, anys ); ASSERT_EQ( val1, A1( 3.5 ) ); ASSERT_EQ( val2, A2( 1.5, 0.5 ) ); ASSERT_EQ( val3, A3( 1.5, 0.5, 0.5 ) );
	multi.get_states( 20, 0, anys ); ASSERT_EQ( val1, A1( 0.5 ) ); ASSERT_EQ( val2, A2( 2.5, 0.5 ) ); ASSERT_EQ( val3, A3( 1.5, 0.5, 0.5 ) );
	multi.get_states( 21, 0, anys ); ASSERT_EQ( val1, A1( 1.5 ) ); ASSERT_EQ( val2, A2( 2.5, 0.5 ) ); ASSERT_EQ( val3, A3( 1.5, 0.5, 0.5 ) );
	multi.get_states( 22, 0, anys ); ASSERT_EQ( val1, A1( 2.5 ) ); ASSERT_EQ( val2, A2( 2.5, 0.5 ) ); ASSERT_EQ( val3, A3( 1.5, 0.5, 0.5 ) );
	multi.get_states( 23, 0, anys ); ASSERT_EQ( val1, A1( 3.5 ) ); ASSERT_EQ( val2, A2( 2.5, 0.5 ) ); ASSERT_EQ( val3, A3( 1.5, 0.5, 0.5 ) );

	#ifndef NDEBUG
	#ifndef CXX14
	ASSERT_DEATH( multi.get_states( 24, 0, anys ), ".*" );
	#endif
	#endif

	multi.get_states(  0 * (1<<6), 1, anys ); ASSERT_EQ( val1, A1( 0.25 ) ); ASSERT_EQ( val2, A2( 0.25, 0.25 ) ); ASSERT_EQ( val3, A3( 0.25, 0.25, 0.25 ) );
	multi.get_states(  1 * (1<<6), 1, anys ); ASSERT_EQ( val1, A1( 1.25 ) ); ASSERT_EQ( val2, A2( 0.25, 0.25 ) ); ASSERT_EQ( val3, A3( 0.25, 0.25, 0.25 ) );
	multi.get_states(  2 * (1<<6), 1, anys ); ASSERT_EQ( val1, A1( 2.25 ) ); ASSERT_EQ( val2, A2( 0.25, 0.25 ) ); ASSERT_EQ( val3, A3( 0.25, 0.25, 0.25 ) );
	multi.get_states(  3 * (1<<6), 1, anys ); ASSERT_EQ( val1, A1( 3.25 ) ); ASSERT_EQ( val2, A2( 0.25, 0.25 ) ); ASSERT_EQ( val3, A3( 0.25, 0.25, 0.25 ) );
	multi.get_states(  4 * (1<<6), 1, anys ); ASSERT_EQ( val1, A1( 0.25 ) ); ASSERT_EQ( val2, A2( 1.25, 0.25 ) ); ASSERT_EQ( val3, A3( 0.25, 0.25, 0.25 ) );
	multi.get_states(  5 * (1<<6), 1, anys ); ASSERT_EQ( val1, A1( 1.25 ) ); ASSERT_EQ( val2, A2( 1.25, 0.25 ) ); ASSERT_EQ( val3, A3( 0.25, 0.25, 0.25 ) );
	multi.get_states(  6 * (1<<6), 1, anys ); ASSERT_EQ( val1, A1( 2.25 ) ); ASSERT_EQ( val2, A2( 1.25, 0.25 ) ); ASSERT_EQ( val3, A3( 0.25, 0.25, 0.25 ) );
	multi.get_states(  7 * (1<<6), 1, anys ); ASSERT_EQ( val1, A1( 3.25 ) ); ASSERT_EQ( val2, A2( 1.25, 0.25 ) ); ASSERT_EQ( val3, A3( 0.25, 0.25, 0.25 ) );
	multi.get_states(  8 * (1<<6), 1, anys ); ASSERT_EQ( val1, A1( 0.25 ) ); ASSERT_EQ( val2, A2( 2.25, 0.25 ) ); ASSERT_EQ( val3, A3( 0.25, 0.25, 0.25 ) );
	multi.get_states(  9 * (1<<6), 1, anys ); ASSERT_EQ( val1, A1( 1.25 ) ); ASSERT_EQ( val2, A2( 2.25, 0.25 ) ); ASSERT_EQ( val3, A3( 0.25, 0.25, 0.25 ) );
	multi.get_states( 10 * (1<<6), 1, anys ); ASSERT_EQ( val1, A1( 2.25 ) ); ASSERT_EQ( val2, A2( 2.25, 0.25 ) ); ASSERT_EQ( val3, A3( 0.25, 0.25, 0.25 ) );
	multi.get_states( 11 * (1<<6), 1, anys ); ASSERT_EQ( val1, A1( 3.25 ) ); ASSERT_EQ( val2, A2( 2.25, 0.25 ) ); ASSERT_EQ( val3, A3( 0.25, 0.25, 0.25 ) );
	multi.get_states( 12 * (1<<6), 1, anys ); ASSERT_EQ( val1, A1( 0.25 ) ); ASSERT_EQ( val2, A2( 0.25, 0.25 ) ); ASSERT_EQ( val3, A3( 1.25, 0.25, 0.25 ) );
	multi.get_states( 13 * (1<<6), 1, anys ); ASSERT_EQ( val1, A1( 1.25 ) ); ASSERT_EQ( val2, A2( 0.25, 0.25 ) ); ASSERT_EQ( val3, A3( 1.25, 0.25, 0.25 ) );
	multi.get_states( 14 * (1<<6), 1, anys ); ASSERT_EQ( val1, A1( 2.25 ) ); ASSERT_EQ( val2, A2( 0.25, 0.25 ) ); ASSERT_EQ( val3, A3( 1.25, 0.25, 0.25 ) );
	multi.get_states( 15 * (1<<6), 1, anys ); ASSERT_EQ( val1, A1( 3.25 ) ); ASSERT_EQ( val2, A2( 0.25, 0.25 ) ); ASSERT_EQ( val3, A3( 1.25, 0.25, 0.25 ) );
	multi.get_states( 16 * (1<<6), 1, anys ); ASSERT_EQ( val1, A1( 0.25 ) ); ASSERT_EQ( val2, A2( 1.25, 0.25 ) ); ASSERT_EQ( val3, A3( 1.25, 0.25, 0.25 ) );
	multi.get_states( 17 * (1<<6), 1, anys ); ASSERT_EQ( val1, A1( 1.25 ) ); ASSERT_EQ( val2, A2( 1.25, 0.25 ) ); ASSERT_EQ( val3, A3( 1.25, 0.25, 0.25 ) );
	multi.get_states( 18 * (1<<6), 1, anys ); ASSERT_EQ( val1, A1( 2.25 ) ); ASSERT_EQ( val2, A2( 1.25, 0.25 ) ); ASSERT_EQ( val3, A3( 1.25, 0.25, 0.25 ) );
	multi.get_states( 19 * (1<<6), 1, anys ); ASSERT_EQ( val1, A1( 3.25 ) ); ASSERT_EQ( val2, A2( 1.25, 0.25 ) ); ASSERT_EQ( val3, A3( 1.25, 0.25, 0.25 ) );
	multi.get_states( 20 * (1<<6), 1, anys ); ASSERT_EQ( val1, A1( 0.25 ) ); ASSERT_EQ( val2, A2( 2.25, 0.25 ) ); ASSERT_EQ( val3, A3( 1.25, 0.25, 0.25 ) );
	multi.get_states( 21 * (1<<6), 1, anys ); ASSERT_EQ( val1, A1( 1.25 ) ); ASSERT_EQ( val2, A2( 2.25, 0.25 ) ); ASSERT_EQ( val3, A3( 1.25, 0.25, 0.25 ) );
	multi.get_states( 22 * (1<<6), 1, anys ); ASSERT_EQ( val1, A1( 2.25 ) ); ASSERT_EQ( val2, A2( 2.25, 0.25 ) ); ASSERT_EQ( val3, A3( 1.25, 0.25, 0.25 ) );
	multi.get_states( 23 * (1<<6), 1, anys ); ASSERT_EQ( val1, A1( 3.25 ) ); ASSERT_EQ( val2, A2( 2.25, 0.25 ) ); ASSERT_EQ( val3, A3( 1.25, 0.25, 0.25 ) );


 }

TEST( MultiNest, get_index_ncell_handling ){
	std::ostream_iterator<uint64_t> out_it (std::cout,", ");
	MultiNest<uint64_t,uint64_t>::Indices indices;

	MultiNest<uint64_t,uint64_t> multi;
	multi.add_nest( make_shared<NEST<1> >(4) );
	multi.add_nest( make_shared<NEST<2> >(3) );
	multi.add_nest( make_shared<NEST<3> >(2) );
	ASSERT_EQ(multi.size(0),24);

	typedef util::SimpleArray<1,double> A1; A1 val1;
	typedef util::SimpleArray<2,double> A2; A2 val2;
	typedef util::SimpleArray<3,double> A3; A3 val3;
	std::vector<boost::any> anys(3);
	anys[0] = &val1;
	anys[1] = &val2;
	anys[2] = &val3;

	for(int resl = 0; resl < 3; ++resl){
		for(int index = 0; index < multi.size(resl); ++index){
			multi.get_states( index, resl, anys );
			// cout << "VALS: (" << val1 << ") (" << val2 << ") (" << val3 << ")"<<endl;
			uint64_t test_index = multi.virtual_get_index( &anys, resl );
			ASSERT_EQ( index, test_index );
		}
	}



 }

}}}
