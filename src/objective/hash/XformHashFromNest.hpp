#ifndef INCLUDED_objective_hash_XformHashFromNest_HH
#define INCLUDED_objective_hash_XformHashFromNest_HH

#include "util/SimpleArray.hpp"

#include <boost/type_traits/is_same.hpp>
#include <boost/static_assert.hpp>

namespace scheme { namespace objective { namespace hash {

template< class _Nest >
struct XformHashFromNest {

	template< class _Xform >
	struct apply {

		typedef boost::is_same< typename _Nest::Value, _Xform > XformCheck;
		BOOST_STATIC_ASSERT_MSG( XformCheck::value, "Nest::Value must be same as Xform" );

		typedef _Nest Nest;
		typedef _Xform Xform;
		typedef typename _Nest::Float Float;
		typedef typename _Nest::Index Index;
		typedef util::SimpleArray<3,Float> F3;
		typedef util::SimpleArray<3,int>   I3;		

		Nest nest_;

		apply(){}

		apply( Float cart_resl, Float ang_resl_deg, Float cart_bound=512.0 ){
			init( cart_resl, ang_resl_deg, cart_bound );
		}

		void
		init( Float cart_resl, Float ang_resl_deg, Float cart_bound ) {
			F3 lb(-cart_bound), ub(cart_bound);
			I3 bs( std::ceil(2.0*cart_bound/cart_resl) );
			nest_.init( ang_resl_deg, lb, ub, bs );
		}

		Index get_key( Xform const & x, int resl = 0 ) const {
			return nest_.get_index( x, resl );
		}

		bool is_valid( Index key, int resl = 0 ) const {
			return nest_.check_state( key, resl );
		}

		Xform get_center( Index key, int resl = 0 ) const {
			Xform val;
			BOOST_VERIFY( nest_.get_state( key, resl, val ) );
			return val;
		}

		Index approx_size( int resl = 0 ) const {
			return nest_.size( resl );
		}


		static std::string name(){ return "XformHashFromNest< "+_Nest::pmap_name() + " >" ; }


	};

};

}}}

#endif
