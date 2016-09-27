#include <gtest/gtest.h>

#include "riflib/objective/hash/XformHashFromNest.hpp"

#include "riflib/nest/pmap/OriTransMap.hpp"
#include "riflib/nest/NEST.hpp"
#include "riflib/objective/hash/XformMap.hpp"


namespace scheme { namespace objective { namespace hash { namespace test_XHFN { 

TEST( XformHashFromNest, static_checks ){

	typedef Eigen::Transform< float, 3, Eigen::AffineCompact > Xform;
	typedef Eigen::Transform< double, 3, Eigen::AffineCompact > Xform2;	
	typedef ::scheme::nest::NEST< 6, Xform, nest::pmap::OriTransMap > Nest;

	// this must fail to compile because Xform2 != Xform
	// typedef XformMap< Xform2, float, XformHashFromNest< Nest >::apply > XMap;
	// XMap xmap;

	typedef XformMap< Xform, float, XformHashFromNest< Nest >::apply > XMap;
	XMap xmap;
}

TEST( XformHashFromNest, construction ){

	typedef Eigen::Transform< float, 3, Eigen::AffineCompact > Xform;
	typedef ::scheme::nest::NEST< 6, Xform, nest::pmap::OriTransMap > Nest;
	typedef XformMap< Xform, float, XformHashFromNest< Nest >::apply > XMap;

	XMap xmap( 1.0, 10.0 );

}

TEST( XformHashFromNest, basic_operation ){

	typedef Eigen::Transform< float, 3, Eigen::AffineCompact > Xform;
	typedef ::scheme::nest::NEST< 6, Xform, nest::pmap::OriTransMap > Nest;
	typedef XformMap< Xform, float, XformHashFromNest< Nest >::apply > XMap;

	XMap xmap( 1.0, 10.0, 100.0 );
	std::cout << xmap.hasher_.nest_.ori_map_.nside_ << std::endl;
	std::cout << xmap.hasher_.nest_.trans_map_.lower_bound() << std::endl;	
	std::cout << xmap.hasher_.nest_.trans_map_.upper_bound() << std::endl;	
	std::cout << xmap.hasher_.nest_.trans_map_.cell_sizes() << std::endl;			

	{
		Xform x( Xform::Identity() );
		// ASSERT_EQ( 9402275379703, xmap.hasher_.get_key( x ) ); // not a real test
	}

	int nfail = 0;

	int const resl  = 0;
	typename Nest::Index end = std::min( (uint64_t)xmap.hasher_.approx_size(resl), (uint64_t)10000 );
	for( typename Nest::Index i = 0; i < end; ++i ){
		// if( ! xmap.hasher_.is_valid( i, r ) ) continue;
		// Xform x = xmap.hasher_.get_center( i, r );
		// ASSERT_EQ( xmap.hasher_.get_key( x, r ), i );
		Xform x;
		if( !xmap.hasher_.nest_.get_state( i, resl, x ) ) continue;
		nfail += xmap.hasher_.nest_.get_index( x, resl ) != i;

	}

	ASSERT_LE( nfail*1./end, 0.05 );

}

}}}}


