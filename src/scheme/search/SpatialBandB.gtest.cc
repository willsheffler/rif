#include <gtest/gtest.h>

#include <scheme/search/SpatialBandB.hh>

#include <scheme/actor/BackboneActor.hh>
#include <scheme/kinematics/Scene.hh>
#include <scheme/objective/ObjectiveFunction.hh>

#include <Eigen/Geometry>

#include <boost/mpl/vector.hpp>

namespace scheme { namespace search { namespace spbbtest {

using std::cout;
using std::endl;

typedef float Float;
typedef Eigen::Transform<Float,3,Eigen::AffineCompact> Xform;
typedef actor::BackboneActor<Xform> BBActor;
typedef Eigen::Matrix<Float,3,1> Vec;


// why is this here?
struct TestScoreBBactor {
	typedef double Result;
	typedef BBActor Interaction;
	static std::string name(){ return "TestScoreBBactor"; }
	Xform target_inv_;
	Float radius_, lever_;
	TestScoreBBactor() : target_inv_(Xform::Identity()), radius_(1.0), lever_(10.0) {}
	TestScoreBBactor( Xform target, Float rad, Float lever ) : target_inv_(target.inverse()), radius_(rad), lever_(lever) {}
	template<class Config>
	Result operator()(BBActor const & a, Config const& c) const {
		// for( int i = 0; i < 12; ++i ){
			// cout << i << " " << target_inv_.data()[i] << endl;
			// cout << i << " " << a.position().data()[i] << endl;
		// }
		// cout << a.position().rotation().trace() << " " << target_inv_.rotation().trace() << endl;
		// cout << target_inv_.inverse().translation().transpose() << endl;
		// cout << a.position().translation().transpose() << endl;
		Xform xrel = target_inv_ * a.position();
		// cout << xrel.translation().transpose() << endl;
		// cout << (xrel.rotation().trace()-1.0)/2.0 << std::endl;
		Float err2 = xrel.translation().squaredNorm();
		Float trace = xrel.data()[0] + xrel.data()[4] +  xrel.data()[8];
		Float roterr = ( 1.0 - ( trace - 1.0 ) / 2.0 ) * lever_;
		// cout << roterr << " " << ( ( trace - 1.0 ) / 2.0 ) << " " << acos(( ( trace - 1.0 ) / 2.0 ) )*180.0/3.14159 << endl;
		roterr = roterr < 0.0 ? 9e9 : roterr;
		Float toterr = err2+roterr*roterr;
		return (toterr < radius_*radius_ ) ? -1.0 : 0.0 ;
	}
};


TEST( SpatialBandB, test_test_utils ){
	typedef kinematics::Scene< boost::mpl::vector< BBActor > , Xform > Scene;
	typedef objective::ObjectiveFunction< boost::mpl::vector< TestScoreBBactor >, int > Objective;

	Objective obj;
	obj.get_objective<TestScoreBBactor>().lever_ = 10.0;
	obj.get_objective<TestScoreBBactor>().radius_ = 1.0;

	{
		Xform x( Xform::Identity() );
		Scene s;
		s.add_body();
		BBActor actor;
		s.add_actor( 0, actor );
		ASSERT_FLOAT_EQ( -1.0, obj(s).sum() );
		{
			x.translation() = Vec(0.99,0,0);
			s.set_position( 0, x );
			ASSERT_FLOAT_EQ( -1.0, obj(s).sum() );
		}
		{
			x.translation() = Vec(1.01,0,0);
			s.set_position( 0, x );
			ASSERT_FLOAT_EQ(  0.0, obj(s).sum() );
		}
		{
			x = Xform::Identity();
			x.rotate( Eigen::AngleAxis<Float>( 25.84/180.0*3.14159, Vec(1,0,0) ) );
			s.set_position( 0, x );
			ASSERT_FLOAT_EQ( -1.0, obj(s).sum() );
		}
		{
			x = Xform::Identity();
			x.rotate( Eigen::AngleAxis<Float>( 25.85/180.0*3.14159, Vec(1,0,0) ) );
			s.set_position( 0, x );
			ASSERT_FLOAT_EQ(  0.0, obj(s).sum() );
		}

		Xform target( Xform::Identity() );
		target.translation() = Vec(1,2,3);
		obj.get_objective<TestScoreBBactor>().target_inv_ = target.inverse();
		{
			x = Xform::Identity();
			x.translation() = Vec(1,2,3);
			s.set_position( 0, x );
			ASSERT_FLOAT_EQ( -1.0, obj(s).sum() );
		}
		{
			x = Xform::Identity();
			x.translation() = Vec(3,2,1);
			s.set_position( 0, x );
			ASSERT_FLOAT_EQ(  0.0, obj(s).sum() );
		}
		{
			x = Xform::Identity();
			x.translation() = Vec(1,2,3.0);
			x.rotate( Eigen::AngleAxis<Float>( 25.84/180.0*3.14159, Vec(1,0,0) ) );
			s.set_position( 0, x );
			ASSERT_FLOAT_EQ( -1.0, obj(s).sum() );
		}
		{
			x = Xform::Identity();
			x.translation() = Vec(1,2,3);
			x.rotate( Eigen::AngleAxis<Float>( 25.85/180.0*3.14159, Vec(1,0,0) ) );
			s.set_position( 0, x );
			ASSERT_FLOAT_EQ(  0.0, obj(s).sum() );
		}
		{
			obj.get_objective<TestScoreBBactor>().lever_ = 60;
			x = Xform::Identity();
			x.translation() = Vec(1,2,3);
			x.rotate( Eigen::AngleAxis<Float>( 25.85/180.0*3.14159, Vec(1,0,0) ) );
			s.set_position( 0, x );
			ASSERT_FLOAT_EQ(  0.0, obj(s).sum() );
		}


	}


	// BoundingObjectiveFunction<


}

}}}
