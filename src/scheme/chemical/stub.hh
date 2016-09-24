#ifndef INCLUDED_chemical_stub_HH
#define INCLUDED_chemical_stub_HH

#include <Eigen/Geometry>

namespace scheme {
namespace chemical {


	template<class Xform,class Point>
	Xform make_stub(Point n, Point ca, Point c){
		Xform out;
		Point cen = ca;
		Point e1( n - ca );
		e1.normalize();
		Point e3( e1.cross(c-ca) );
		e3.normalize();
		Point e2( e3.cross(e1) );
		out.linear().col(0) = e1;
		out.linear().col(1) = e2;
		out.linear().col(2) = e3;
		out.translation() = cen;
		return out;
	}


}
}

#endif
