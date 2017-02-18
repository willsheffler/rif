#ifndef INCLUDED_scheme_nest_maps_UnitMap_HH
#define INCLUDED_scheme_nest_maps_UnitMap_HH

#include "util/template_loop.hpp"
#include "util/SimpleArray.hpp"
#include <boost/function.hpp>
#include <boost/bind.hpp>
#include <boost/static_assert.hpp>
#include <iostream>
#include <vector>

namespace scheme {
namespace nest {
namespace pmap {

	///@brief Parameter to Value Map Policy Class
	///@detail just copies the [0.0,1.0] hypercube coords to Value
	///        the first dimension will incremented by cell_index
	///@tparam DIM the dimension number of the input parameter space
	///@tparam Value the output value type, default SimpleArray
	///@tparam Index index type, default size_t
	///@tparam Float float type, default double
	template<
		int DIM,
		class Value=util::SimpleArray<DIM,double>,
		class Index=uint64_t,
		class Float=double
	>
	struct UnitMap {
		static int const DIMENSION = DIM;
		typedef Value ValueType ;
		typedef Float FloatType ;		
		typedef Index IndexType ;		
		typedef util::SimpleArray<DIM,Index> Indices;
		typedef util::SimpleArray<DIM,Float> Params;

		BOOST_STATIC_ASSERT_MSG(DIM>0,"UnitMap DIM must be > 0");
		///@brief constructor
		UnitMap(Index num_cells=1) : num_cells_(num_cells) {}
		///@brief sets value to parameters without change
		///@return false iff invalid parameters
		bool params_to_value(
			Params const & params,
			Index cell_index,
			Index /*resl*/,
			Value & value
		) const {
			for(int i = 0; i < DIM; ++i) assert( 0.0 <= params[i] );
			assert( params[0] <= (Float)num_cells_ );
			for(int i = 1; i < DIM; ++i) assert( params[i] <= 1.0 );
            for(int i = 0; i < DIM; ++i) value[i] = params[i];
			value[0] += (Float)cell_index;
			return true;
		}
		///@brief sets params/cell_index from value
		///@note necessary for value lookup and neighbor lookup
		bool value_to_params(
			Value const & value,
			Index resl,
			Params & params,
			Index & cell_index
		) const {
			///@note neighbor lookups require out of bounds mappings to be valid
			value_to_params_for_cell(value,resl,params,0);
			cell_index = (Index)value[0];
			params[0] -= (Float)cell_index;
			for(int i = 0; i < DIM; ++i) assert( 0.0 <= params[i] );
			assert( params[0] <= (Float)num_cells_ );
			for(int i = 1; i < DIM; ++i) assert( params[i] <= 1.0 );
			return true;
		}
		///@brief get params repr of Value wrt cell cell_index
		///@note necessary only for neighbor lookup		
		void value_to_params_for_cell(
			Value const & value,
			Index /*resl*/,
			Params & params,
			Index cell_index
		) const {
			for(int i = 0; i < DIM; ++i) params[i] = value[i];
			params[0] -= (Float)cell_index;
		}
		///@brief return the cell_index of neighboring cells within delta of value
		///@note delta parameter is in "Parameter Space"
		template<class OutIter>
		void get_neighboring_cells(
			Value const & value,
			Index /*resl*/,
			Float param_delta,
			OutIter out
		) const {
			// Float param_delta = 1.0 / (Float)(1<<resl);
			// this BIG thing is to ensure rounding goes down
			int const BIG = 12345678;
			int lb = std::max(                0, static_cast<int>( value[0]-param_delta + BIG ) - BIG );
			int ub = std::min((int)num_cells_-1, static_cast<int>( value[0]+param_delta + BIG ) - BIG );
			// std::cout << "lb " << lb << " ub "  << ub << std::endl;
			// assert(lb<=ub);
			for(int i = lb; i <= ub; ++i) *(out++) = i;
		}
		///@brief aka covering radius max distance from bin center to any value within bin
		Float bin_circumradius(Index resl) const { return 0.5/(Float)(1<<resl) * sqrt(DIM); }
		///@brief maximum distance from the bin center which must be within the bin
		Float bin_inradius(Index resl) const { return 1.5/(Float)(1<<resl); }
		///@brief cell size
		Index num_cells() const { return num_cells_; }
		virtual ~UnitMap(){}
	 private:
	 	///@brief number of cells
		Index num_cells_;
	};
	


}
}
}

#endif
