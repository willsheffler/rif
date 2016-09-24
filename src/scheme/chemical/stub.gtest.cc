#include <gtest/gtest.h>

#include "scheme/chemical/stub.hh"

#include <iostream>

namespace scheme {
namespace chemical {
namespace test234 {

using namespace Eigen;
typedef ::Eigen::Transform<float,3,Eigen::AffineCompact> Xform;


TEST( stub, basic_alignment ){
	Vector3f a1(1,0,0), b1(0,0,0), c1(0,1,0);
	Vector3f a2(1,2,1), b2(1,1,1), c2(1,1,2);


	Xform x1 = make_stub<Xform>(a1,b1,c1);
	std::cout << x1.linear() << std::endl;
	std::cout << x1.translation().transpose() << std::endl;
	std::cout << std::endl;
	Xform x2 = make_stub<Xform>(a2,b2,c2);
	std::cout << x2.linear() << std::endl;
	std::cout << x2.translation().transpose() << std::endl;

	Xform x1to2 = x2*x1.inverse();
	ASSERT_LT( (x1to2*a1 - a2).norm() , 0.0001 );
	ASSERT_LT( (x1to2*b1 - b2).norm() , 0.0001 );
	ASSERT_LT( (x1to2*c1 - c2).norm() , 0.0001 );

	Xform x2to1 = x1*x2.inverse();
	ASSERT_LT( (x2to1*a2 - a1).norm() , 0.0001 );
	ASSERT_LT( (x2to1*b2 - b1).norm() , 0.0001 );
	ASSERT_LT( (x2to1*c2 - c1).norm() , 0.0001 );

}


}
}
}

