#include <gtest/gtest.h>

#include "numeric/util.hpp"

#include "nest/pmap/Rotation1DMap.hpp"
#include "nest/NEST.hpp"

#include <Eigen/Dense>
#include <random>

#include <boost/lexical_cast.hpp>

#include <fstream>

namespace scheme { namespace nest { namespace pmap { namespace test {

using namespace Eigen;
using std::cout;
using std::endl;

Vector3d get_axis( Matrix3d const & m ){
	AngleAxisd aa( m );
	if( aa.axis().dot( Vector3d(1,1,1) ) < 0 ) return -aa.axis();
	else                                       return  aa.axis();
}
double get_angle( Matrix3d const & m ){
	AngleAxisd aa( m );
	if( aa.axis().dot( Vector3d(1,1,1) ) < 0 ) return -aa.angle();
	else                                       return  aa.angle();
}

TEST( Rotation1DMap , flip_test ){

	{
		NEST<1,Matrix3d,Rotation1DMap> nest;
		Vector3d axis(1,0,0), flip(0,1,0);
		nest.set_axis( axis );
		nest.set_flip_axis( flip );
		{ Matrix3d m = nest.set_and_get(0,0); ASSERT_DOUBLE_EQ(  0.0     , get_angle(m) );                                               }
		{ Matrix3d m = nest.set_and_get(1,0); ASSERT_DOUBLE_EQ(  M_PI    , get_angle(m) ); ASSERT_GE( get_axis(m).dot( flip ), 0.9999 ); }
		{ Matrix3d m = nest.set_and_get(0,1); ASSERT_DOUBLE_EQ( -M_PI/2.0, get_angle(m) ); ASSERT_GE( get_axis(m).dot( axis ), 0.9999 ); }
		{ Matrix3d m = nest.set_and_get(1,1); ASSERT_DOUBLE_EQ(  M_PI/2.0, get_angle(m) ); ASSERT_GE( get_axis(m).dot( axis ), 0.9999 ); }
		{
			Matrix3d m = nest.set_and_get(2,1);
			ASSERT_DOUBLE_EQ(  M_PI    , get_angle(m) );
			ASSERT_GE( get_axis(m).dot( Vector3d(0,sqrt(2)/2,-sqrt(2)/2) ), 0.9999 );
		}
		{
			Matrix3d m = nest.set_and_get(3,1);
			ASSERT_DOUBLE_EQ(  M_PI, get_angle(m) );
			ASSERT_GE( get_axis(m).dot( Vector3d(0,sqrt(2)/2,sqrt(2)/2) ), 0.9999 );
		}
		#ifndef NDEBUG
		#ifndef CXX14
		ASSERT_DEATH( nest.set_and_get(4,1), ".*" );
		#endif
		#endif
	}

}

TEST( Rotation1DMap , basic_test ){

	{
		NEST<1,Matrix3d,Rotation1DMap> nest;
		Vector3d axis(1,0,0); nest.set_axis( axis );
		{ Matrix3d m = nest.set_and_get(0,0); ASSERT_DOUBLE_EQ(  0.0     , get_angle(m) );                                               }
		{ Matrix3d m = nest.set_and_get(0,1); ASSERT_DOUBLE_EQ( -M_PI/2.0, get_angle(m) ); ASSERT_GE( get_axis(m).dot( axis ), 0.9999 ); }
		{ Matrix3d m = nest.set_and_get(1,1); ASSERT_DOUBLE_EQ(  M_PI/2.0, get_angle(m) ); ASSERT_GE( get_axis(m).dot( axis ), 0.9999 ); }
	}
	{
		NEST<1,Matrix3d,Rotation1DMap> nest;
		Vector3d axis(1,1,1); nest.set_axis( axis );
		{ Matrix3d m = nest.set_and_get(0,0); ASSERT_DOUBLE_EQ(  0.0     , get_angle(m) );                                               }
		{ Matrix3d m = nest.set_and_get(0,1); ASSERT_DOUBLE_EQ( -M_PI/2.0, get_angle(m) ); ASSERT_GE( get_axis(m).dot( axis ), 0.9999 ); }
		{ Matrix3d m = nest.set_and_get(1,1); ASSERT_DOUBLE_EQ(  M_PI/2.0, get_angle(m) ); ASSERT_GE( get_axis(m).dot( axis ), 0.9999 ); }
	}
	{
		NEST<1,Matrix3d,Rotation1DMap> nest(0,M_PI,1);
		Vector3d axis(1,1,1); nest.set_axis( axis );
		{ Matrix3d m = nest.set_and_get(0,0); ASSERT_DOUBLE_EQ(    M_PI/2.0, get_angle(m) ); ASSERT_GE( get_axis(m).dot( axis ), 0.9999 ); }
		{ Matrix3d m = nest.set_and_get(0,1); ASSERT_DOUBLE_EQ(    M_PI/4.0, get_angle(m) ); ASSERT_GE( get_axis(m).dot( axis ), 0.9999 ); }
		{ Matrix3d m = nest.set_and_get(1,1); ASSERT_DOUBLE_EQ(  3*M_PI/4.0, get_angle(m) ); ASSERT_GE( get_axis(m).dot( axis ), 0.9999 ); }
	}


}


TEST( Rotation1DMap , lookup ){
	int NITER = 1000;
	#ifdef NDEBUG
	NITER *= 30;
	#endif

	NEST<1,Matrix3d,Rotation1DMap> nest;
	// cout << nest.get_index( Matrix3d::Identity(), 0 ) << endl;

	NestBase<> *nestp = new NEST<1,Matrix3d,Rotation1DMap>();
	Matrix3d m = Matrix3d::Identity();
	ASSERT_EQ( 0, nestp->virtual_get_index( &m , 0 ) );

	std::mt19937 rng(0);
	std::normal_distribution<> gauss;
	std::uniform_real_distribution<> uniform;
	for(int i = 0; i < NITER; ++i){
		Vector3d axis( gauss(rng), gauss(rng), gauss(rng) );
		axis.normalize();
		nest.set_axis(axis);
		nest.lb = uniform(rng) * 2*M_PI - M_PI;
		nest.ub = nest.lb + uniform(rng) * (2*M_PI-nest.lb-M_PI);
		BOOST_VERIFY( -M_PI <= nest.lb && nest.lb <= M_PI );
		BOOST_VERIFY( -M_PI <= nest.ub && nest.ub <= M_PI );
		BOOST_VERIFY( nest.lb < nest.ub );
		nest.nside = (uint64_t)(3*uniform(rng)+1.0);

		for(int resl = 0; resl < 5; ++resl){
			for(uint64_t index = 0; index < nest.size(resl); ++index){
				Matrix3d m = nest.set_and_get(index,resl);
				ASSERT_EQ( nest.get_index(m,resl), index );
			}
		}


	}
}

TEST( Rotation1DMap , lookup_flip ){
	int NITER = 1000;
	#ifdef NDEBUG
	NITER *= 30;
	#endif

	NEST<1,Matrix3d,Rotation1DMap> nest;
	// cout << nest.get_index( Matrix3d::Identity(), 0 ) << endl;

	NestBase<> *nestp = new NEST<1,Matrix3d,Rotation1DMap>();
	Matrix3d m = Matrix3d::Identity();
	ASSERT_EQ( 0, nestp->virtual_get_index( &m , 0 ) );

	std::mt19937 rng(0);
	std::normal_distribution<> gauss;
	std::uniform_real_distribution<> uniform;
	for(int i = 0; i < NITER; ++i){
		Vector3d axis( gauss(rng), gauss(rng), gauss(rng) );
		axis.normalize();
		nest.set_axis(axis);
		nest.set_flip_axis( axis.cross( Vector3d( gauss(rng), gauss(rng), gauss(rng) ) ) );
		nest.lb = uniform(rng) * 2*M_PI - M_PI;
		nest.ub = nest.lb + uniform(rng) * (2*M_PI-nest.lb-M_PI);
		BOOST_VERIFY( -M_PI <= nest.lb && nest.lb <= M_PI );
		BOOST_VERIFY( -M_PI <= nest.ub && nest.ub <= M_PI );
		BOOST_VERIFY( nest.lb < nest.ub );
		nest.nside = (uint64_t)(3*uniform(rng)+1.0);

		for(int resl = 0; resl < 5; ++resl){
			for(uint64_t index = 0; index < nest.size(resl); ++index){
				Matrix3d m = nest.set_and_get(index,resl);
				ASSERT_EQ( nest.get_index(m,resl), index );
			}
		}


	}
}

}}}}
