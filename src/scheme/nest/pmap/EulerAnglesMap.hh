#ifndef INCLUDED_scheme_nest_maps_EulerAnglesMap_HH
#define INCLUDED_scheme_nest_maps_EulerAnglesMap_HH

#include "scheme/util/SimpleArray.hh"
#include "scheme/numeric/euler_angles.hh"

#include <boost/static_assert.hpp>
// #include <iostream>
#include <vector>

namespace scheme {
namespace nest {
namespace pmap {



template<
	int DIM,
	class Value=util::SimpleArray<DIM,double>,
	class Index=size_t,
	class Float=double
>
struct EulerAnglesMap {
	static int const DIMENSION = DIM;
	typedef Value ValueType ;
	typedef Float FloatType ;		
	typedef Index IndexType ;		
	typedef util::SimpleArray<DIM,Index> Indices;
	typedef util::SimpleArray<DIM,Float> Params;

	BOOST_STATIC_ASSERT_MSG(DIM==3,"EulerAnglesMap DIM must be == 3");
	///@brief constructor
	EulerAnglesMap() {}
	///@brief
	///@return false iff invalid parameters
	bool params_to_value(
		Params const & params,
		Index /*cell_index*/,
		Index /*resl*/,
		Value & value
	) const {
		for(size_t i = 0; i < DIM; ++i) assert( 0.0 <= params[i] );
		// assert( params[0] <= (Float)num_cells_ );
		for(size_t i = 1; i < DIM; ++i) assert( params[i] <= 1.0 );
		if( params[2] > 0.5 ) return false; // convention is 'W' component is >= 0
		Params euler = params*2.0*boost::math::constants::pi<Float>();
		numeric::from_euler_angles(euler,value);
		return true;
	}
	///@brief
	///@note necessary for value lookup and neighbor lookup
	bool value_to_params(
		Value const & value,
		Index resl,
		Params & params,
		Index & cell_index
	) const {
		///@note neighbor lookups require out of bounds mappings to be valid
		cell_index = 0;
		value_to_params_for_cell(value,resl,params,0);
		return true;
	}
	///@brief
	///@note necessary only for neighbor lookup		
	void value_to_params_for_cell(
		Value const & value,
		Index /*resl*/,
		Params & params,
		Index /*cell_index*/
	) const {
		numeric::euler_angles(value,params);
		// std::cout << params << std::endl;
		params = params / 2.0 / boost::math::constants::pi<Float>();
		for(size_t i = 0; i < DIM; ++i) assert( 0.0 <= params[i] );
		for(size_t i = 1; i < DIM; ++i) assert( 1.0 >= params[i] );
	}
	///@brief
	///@note delta parameter is in "Parameter Space"
	template<class OutIter>
	void get_neighboring_cells(
		Value const & value,
		Index /*resl*/,
		Float param_delta,
		OutIter out
	) const {
		// // Float param_delta = 1.0 / (Float)(1<<resl);
		// // this BIG thing is to ensure rounding goes down
		// int const BIG = 12345678;
		// int lb = std::max(                0, static_cast<int>( value[0]-param_delta + BIG ) - BIG );
		// int ub = std::min((int)num_cells_-1, static_cast<int>( value[0]+param_delta + BIG ) - BIG );
		// // std::cout << "lb " << lb << " ub "  << ub << std::endl;
		// // assert(lb<=ub);
		// for(int i = lb; i <= ub; ++i) *(out++) = i;
	}
	///@brief aka covering radius max distance from bin center to any value within bin
	Float bin_circumradius(Index resl) const { return 2.0/(Float)(1<<resl); }
	///@brief maximum distance from the bin center which must be within the bin
	Float bin_inradius(Index resl) const { return 3.0/(Float)(1<<resl); }
	///@brief cell size
	Index num_cells() const { return 1; }
	virtual ~EulerAnglesMap(){}
 private:
};



}
}
}

#endif
