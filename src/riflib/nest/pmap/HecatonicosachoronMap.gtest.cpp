#include <gtest/gtest.h>

#include "riflib/nest/NEST.hpp"
#include "HecatonicosachoronMap.hpp"
#include "riflib/util/str.hpp"
#include "riflib/io/dump_pdb_atom.hpp"

#include <Eigen/Geometry>
// #include <Eigen/Dense>
#include <iostream>
#include <fstream>
#include <vector>
#include <boost/foreach.hpp>
#include <random>
#include <boost/lexical_cast.hpp>

namespace scheme { namespace nest { namespace pmap { namespace hecat_test {

using namespace Eigen;
using std::cout;
using std::endl;
using std::setprecision;
using std::vector;
using std::ofstream;
using std::string;
// using std::bitset;
using std::ostream;


double const r = (1.0+sqrt(5.0))/2.0;

void softassert( bool b, string s ){
	if(!b) cout << "WARNING: " << s << endl;
}

void
dumpatom(
	ostream & out,
	bool nothet,
	size_t iatom,
	string atom,
	string res,
	char chain,
	size_t ires,
	double x, double y, double z,
	double o, double b,
	string elem
){
	assert( atom.size()<5);
	assert( res.size()<4);
	softassert( x<10000 && x > -1000, "x out of bounds" );
	softassert( y<10000 && y > -1000, "y out of bounds" );
	softassert( z<10000 && z > -1000, "z out of bounds" );		
	// cout << "ATOM   1604  C   GLU A 220       5.010  12.933   1.553  1.00 41.10           C" << endl;
	char buf[128];
	snprintf(buf,128,"%s%5lu %4s %3s %c%4lu    %8.3f%8.3f%8.3f%6.2f%6.2f %11s\n",
		nothet?"ATOM  ":"HETATM",
		iatom,
		atom.c_str(),
		res.c_str(),
		chain,
		ires,
		x,y,z
		,o,b,elem.c_str()
	);
	out << buf;
}


ostream & operator<<(ostream & out,Quaterniond const & q){
	char buf[128];
	snprintf(buf,128,"Quaterniond( %8.5f, %8.5f, %8.5f, %8.5f )",q.w(),q.x(),q.y(),q.z());
	out << buf;
	return out;
}

// Matrix4d projmat( Quaterniond q ){
// 	return Vector4d( q.w(), q.x(), q.y(), q.z() ) * RowVector4d( q.w(), q.x(), q.y(), q.z() );
// }

void dumpquats(string fname,vector<Quaterniond> &quats, Matrix4d proj, Quaterniond dir, double window=0.1){
	dir.normalize();
	ofstream out(fname.c_str());
	size_t ia=0, ir=0;
	int count = 0;
	BOOST_FOREACH(Quaterniond q,quats){
		// if( dir.dot(q) < mindot ) continue;
		// if( dir.dot(q) < 0.99 ) continue;
		Vector4d wxyz = proj*Vector4d(q.w(),q.x(),q.y(),q.z());
		if( fabs(wxyz[1]) > window || fabs(wxyz[2]) > window || fabs(wxyz[3]) > window ) continue;
		char chain = 'A';
		if( q.x() >= q.w() && q.x() >= q.x() && q.x() >= q.y() && q.x() >= q.z() ) chain = 'B';
		if( q.y() >= q.w() && q.y() >= q.x() && q.y() >= q.y() && q.y() >= q.z() ) chain = 'C';
		if( q.z() >= q.w() && q.z() >= q.x() && q.z() >= q.y() && q.z() >= q.z() ) chain = 'D';
		dumpatom(out,true,++ia,"QUAT","QUT",chain,++ir,3000*wxyz[1],3000*wxyz[2],3000*wxyz[3],1.0,1.0,"C");
		++count;
	}
	cout << dir << " dumped " << count << endl;
	out.close();
}

template<class M>
void dumppoints( string fname, M const & pts ){
	ofstream out(fname.c_str());
	size_t ia=0, ir=0;
	int count = 0;
	for(int i = 0; i < pts.cols(); ++i){
		// if( dir.dot(q) < mindot ) continue;
		dumpatom(out,true,++ia,"PNT ","PNT",'A',++ir,100.0*pts(0,i),100.0*pts(1,i),100.0*pts(2,i),1.0,1.0,"C");
		++count;
	}
	out.close();
}










// double sign(double d){
// 	if( d == 0) return 0;
// 	if( d > 0 ) return 1.0;
// 	else return -1.0;
// }

template<class T>
double min_cube_side(T const & coords, Matrix<typename T::Scalar,3,1> & cen){
	Matrix<typename T::Scalar,3,1> const mn = coords.rowwise().minCoeff();
	Matrix<typename T::Scalar,3,1> const mx = coords.rowwise().maxCoeff();
	Matrix<typename T::Scalar,3,1> const side = mx-mn;
	cen = (mn+mx)/2.0;
	double const cube = side.maxCoeff();
	return cube;
}
template<class T>
double min_cube_side(T const & coords){
	Matrix<typename T::Scalar,3,1> cen;
	return min_cube_side(coords,cen);
}

// TEST(hecatonicosachoron, sample_dodec_inscribe_cube){
// 	// dodecahedron 
// 	// (±1, ±1, ±1)
// 	// (0, ±1/φ, ±φ)
// 	// (±1/φ, ±φ, 0)
// 	// (±φ, 0, ±1/φ)
// 	Matrix<double,3,20> dodec0;
// 	dodec0 << Vector3d(     +1,   +1,    +1    ).normalized(),
// 	          Vector3d(     +1,   +1,    -1    ).normalized(),
// 	          Vector3d(     +1,   -1,    +1    ).normalized(),
// 	          Vector3d(     +1,   -1,    -1    ).normalized(),
// 	          Vector3d(     -1,   +1,    +1    ).normalized(),
// 	          Vector3d(     -1,   +1,    -1    ).normalized(),
// 	          Vector3d(     -1,   -1,    +1    ).normalized(),
// 	          Vector3d(     -1,   -1,    -1    ).normalized(),
// 	          Vector3d(      0, +1.0/r,     +r ).normalized(),
// 	          Vector3d(      0, +1.0/r,     -r ).normalized(),
// 	          Vector3d(      0, -1.0/r,     +r ).normalized(),
// 	          Vector3d(      0, -1.0/r,     -r ).normalized(),
// 	          Vector3d( +1.0/r,     +r,    0   ).normalized(),
// 	          Vector3d( +1.0/r,     -r,    0   ).normalized(),
// 	          Vector3d( -1.0/r,     +r,    0   ).normalized(),
// 	          Vector3d( -1.0/r,     -r,    0   ).normalized(),
// 	          Vector3d(     +r,    0  , +1.0/r ).normalized(),
// 	          Vector3d(     +r,    0  , -1.0/r ).normalized(),
// 	          Vector3d(     -r,    0  , +1.0/r ).normalized(),
// 	          Vector3d(     -r,    0  , -1.0/r ).normalized();
// 	// dumppoints("dodec0.pdb",dodec0);

// 	std::mt19937 mt((unsigned int)time(0));
// 	std::normal_distribution<> rnorm;
 
// 	 // 0.90919
// 	Quaterniond Q0( 
// 			-0.18683389537307968,
// 			-0.94619545649322123,
// 			0.1868233140277023,
// 			-0.18682693324842686
// 	);
// 	double const rad = 0.001;

// 	Vector3d cen;
// 	double const side0 = min_cube_side(dodec0,cen);
// 	double const area0 = side0*side0*side0;
// 	cout << "min cube side: " << side0 << endl;
// 	BOOST_VERIFY( cen.isApprox(Vector3d(0,0,0)));
// 	double minside = side0;
// 	Matrix<double,3,20> bestdodec;
// 	Quaterniond bestrot;
// 	for(size_t i = 0; i < 100000; ++i){
// 		// Quaterniond qrand = Quaterniond( rnorm(mt), rnorm(mt), rnorm(mt), rnorm(mt) );
// 		Quaterniond qrand( Q0.w()+rnorm(mt)*rad, Q0.x()+rnorm(mt)*rad, Q0.y()+rnorm(mt)*rad, Q0.z()+rnorm(mt)*rad );
// 		qrand.normalize();
// 		if(!i) qrand = Quaterniond(1,0,0,0);
// 		Matrix<double,3,20> rotdodec = qrand.matrix() * dodec0;
// 		double side = min_cube_side(rotdodec,cen);
// 		assert( cen.isApprox(Vector3d(0,0,0)));
// 		if( minside > side ){
// 			minside = side;
// 			bestrot = qrand;
// 			bestdodec = rotdodec;
// 			double const area = side*side*side;
// 			cout << "min cube side: " << side << " vol: " << area << " dodec vol frac: " << area / area0 << " " << i << endl;
// 			cout << bestrot << endl;
// 		}
// 	}
// 	cout << "BEST ROTATION OF DODEC: " << endl; 
// 	cout << std::setprecision(17) << "\t\t\t" << bestrot.w() << ",\n\t\t\t" << bestrot.x() << ",\n\t\t\t" << bestrot.y() << ",\n\t\t\t" << bestrot.z() << endl;
// 	// dodec volume as constructed 2.78516386312, from this:
// 	// print Vec(1,1,1).normalized().distance( Vec(1/r, r, 0).normalized() )**3 * (15+7*sqrt(5)) /4.0


// 	dumppoints("0_dodec_0.pdb",dodec0/dodec0.block(0,0,3,1).norm());
// 	dumppoints("0_dodec_best.pdb",bestdodec/bestdodec.block(0,0,3,1).norm());
// 	// std::exit(0);

// }


vector<Quaterniond> make_half120cell(){
	// hexacosichoron (600-cell), 120 virts
	// (±½,±½,±½,±½),
	// and 8 vertices obtained from
	// (0,0,0,±1)
	// by permuting coordinates. The remaining 96 vertices are obtained by taking even permutations of
	// ½(±φ,±1,±1/φ,0).

	// hecatonicosachoron (120-cell), 600 virts
	// The 600 vertices of the 120-cell include all permutations of:[1]
	// (0, 0, ±2, ±2)
	// (±1, ±1, ±1, ±√5)
	// (±φ-2, ±φ, ±φ, ±φ)
	// (±φ-1, ±φ-1, ±φ-1, ±φ2)
	// and all even permutations of
	// (0, ±φ-2, ±1, ±φ2)
	// (0, ±φ-1, ±φ, ±√5)
	// (±φ-1, ±1, ±φ, ±2)
	// where φ (also called τ) is the golden ratio, (1+√5)/2.
	vector<Quaterniond> half120cell;
	{
		double const r = (1.0+sqrt(5.0))/2.0;
		vector<Quaterniond> c120;
		for(int s1 = 1; s1 > -2; s1-=2){
			// and 8 vertices obtained from (0,0,0,±1)
			c120.push_back( Quaterniond( s1, 0, 0, 0 ) );
			c120.push_back( Quaterniond(  0,s1, 0, 0 ) );
			c120.push_back( Quaterniond(  0, 0,s1, 0 ) );
			c120.push_back( Quaterniond(  0, 0, 0,s1 ) );
		for(int s2 = 1; s2 > -2; s2-=2){
		for(int s3 = 1; s3 > -2; s3-=2){
			// (±½,±½,±½,±½),
			c120.push_back( Quaterniond(  0.5, s1*0.5, s2*0.5, s3*0.5 ) );
			c120.push_back( Quaterniond( -0.5, s1*0.5, s2*0.5, s3*0.5 ) );
			// The remaining 96 vertices are obtained by taking even permutations of
			// ½(±φ,±1,±1/φ,0).
			double const a = s1*r/2.0, b=s2*1/2.0, c=s3/r/2.0, d=0;
			c120.push_back( Quaterniond( a,b,c,d ) ); // 0
			// c120.push_back( Quaterniond( a,b,d,c ) ); // 1
			// c120.push_back( Quaterniond( a,c,b,d ) ); // 1
			c120.push_back( Quaterniond( a,c,d,b ) ); // 2
			c120.push_back( Quaterniond( a,d,b,c ) ); // 2
			// c120.push_back( Quaterniond( a,d,c,b ) ); // 3
			// c120.push_back( Quaterniond( b,a,c,d ) ); // 1
			c120.push_back( Quaterniond( b,a,d,c ) ); // 2
			c120.push_back( Quaterniond( b,c,a,d ) ); // 2
			// c120.push_back( Quaterniond( b,c,d,a ) ); // 3
			// c120.push_back( Quaterniond( b,d,a,c ) ); // 3
			c120.push_back( Quaterniond( b,d,c,a ) ); // 4
			c120.push_back( Quaterniond( c,a,b,d ) ); // 2
			// c120.push_back( Quaterniond( c,a,d,b ) ); // 3
			// c120.push_back( Quaterniond( c,b,a,d ) ); // 3
			c120.push_back( Quaterniond( c,b,d,a ) ); // 4
			c120.push_back( Quaterniond( c,d,a,b ) ); // 4
			// c120.push_back( Quaterniond( c,d,b,a ) ); // 5
			// c120.push_back( Quaterniond( d,a,b,c ) ); // 3
			c120.push_back( Quaterniond( d,a,c,b ) ); // 4
			c120.push_back( Quaterniond( d,b,a,c ) ); // 4
			// c120.push_back( Quaterniond( d,b,c,a ) ); // 5
			// c120.push_back( Quaterniond( d,c,a,b ) ); // 5
			c120.push_back( Quaterniond( d,c,b,a ) ); // 6
		}}}


		/// provieds about 10% decrease to out-of-bounds samples
		Quaterniond alignq( 
			-0.18683389537307968,
			-0.94619545649322123,
			 0.1868233140277023,
	 		-0.18682693324842686
		);
		Matrix4d alignment; alignment.fill(0);
		alignment.block(0,0,3,3) = alignq.matrix();
		alignment(3,3) = 1.0;
		BOOST_VERIFY( Vector4d(0,0,0,1).isApprox( alignment*Vector4d(0,0,0,1) ) );

		BOOST_FOREACH(Quaterniond q,c120){
			q.coeffs() = alignment * q.coeffs();
			if(  q.w()> 0 || 
				(q.w()==0 && q.x()>0) || 
				(q.w()==0 && q.x()==0 && q.y()>0) || 
				(q.w()==0 && q.x()==0 && q.y()==0 && q.z()>0)
			){
				half120cell.push_back( q );
			// } else {
				// BOOST_VERIFY( false );
			}
		}

		if( half120cell.size() != 60 ){
			cout << "problem making half-120-cell " << half120cell.size() << endl;
			// std::exit(-1);
		}
	}

	return half120cell;
}


Array<uint8_t,60,12>
make_half120cell_neighbors(
	vector<Quaterniond> const & half120cell
){
	Array<uint8_t,60,12> nbrs;
	for(size_t i = 0; i < 60; ++i){
		Quaterniond const & p = half120cell[i];
		size_t nnbr = 0;
		for(size_t j = 0; j < 60; ++j){
			if( i==j ) continue;
			Quaterniond const & q = half120cell[j];
			double ang = acos( 2.0*q.dot(p)*q.dot(p) - 1.0 )*180.0/M_PI;
			if( ang < 73.0 ){
				BOOST_VERIFY( nnbr < 12 );
				nbrs(i,nnbr++) = (uint8_t)j;
			}
		}
		BOOST_VERIFY( nnbr == 12 );
	}
	// for(size_t i = 0; i < 60; ++i){
	// 	cout << "nbrs " << i;
	// 	for(size_t j = 0; j < 12; ++j){
	// 		cout << " " << (int)nbrs(i,j);
	// 	}
	// 	cout << endl;
	// }
	return nbrs;
}
vector<Quaterniond> make_half120cell();


// TEST(hecatonicosachoron, dump_static_data){
// 	vector<Quaterniond> half120cell = make_half120cell();
// 	Array<uint8_t,60,12> nbrs = make_half120cell_neighbors(half120cell);

// 	cout << "\t\tstatic T const h120[240] = {" << endl;
// 	BOOST_FOREACH(Quaterniond q,half120cell){
// 		cout << "\t\t\t" << setprecision(17)  << q.x() << "," << q.y() << "," << q.z() << "," << q.w() << "," << endl;
// 	}
// 	cout << "\t\t}" << endl;

// 	cout << "\t\tstatic T const h120inv[240] = {" << endl;
// 	BOOST_FOREACH(Quaterniond q,half120cell){
// 		q = q.inverse();
// 		cout << "\t\t\t" << setprecision(17)  << q.x() << "," << q.y() << "," << q.z() << "," << q.w() << "," << endl;
// 	}
// 	cout << "\t\t}" << endl;

// 	cout << "\t\tstatic uint8_t const h120_nbrs[60,12] = {" << endl;
// 	for(int i = 0; i < 60; ++i){
// 		cout << "\t\t\t";
// 		for(int j = 0; j < 12; ++j){		
// 			cout << static_cast<int>(nbrs(i,j)) << ",";
// 		}
// 		cout << endl;
// 	}
// 	cout << "\t\t}" << endl;

// 	cout << "\t\tstatic T const cellfaces[36] = {" << endl;
// 	for(int i = 0; i < 12; ++i){
// 		Vector3d nbr0 = half120cell[nbrs(0,i)].coeffs().block(0,0,3,1);
// 		nbr0 /= nbr0.norm();
// 		cout << "\t\t\t" << setprecision(17)  << nbr0[0] << "," << nbr0[1] << "," << nbr0[2] << "," << endl;
// 	}
// 	cout << "\t\t}" << endl;


// }

// TEST(hecatonicosachoron, cell_alignment){
// 	/////////////////////////////////////////////////////////////////////////////////////
// 	/////////////////////////////////////////////////////////////////////////////////////
// 	/////////////////////////////////////////////////////////////////////////////////////		
// 	///////////// this doesn't make sense, must use vertices not centers! ///////////////
// 	/////////////////////////////////////////////////////////////////////////////////////
// 	/////////////////////////////////////////////////////////////////////////////////////		
// 	vector<Quaterniond> half120cell = make_half120cell();
// 	Array<uint8_t,60,12> nbrs = make_half120cell_neighbors(half120cell);

// 	std::mt19937 mt((unsigned int)time(0));
// 	std::normal_distribution<> rnorm;
 
// 	Matrix<double,4,13> dodec0;
// 	for(size_t inbr = 0; inbr < 12; ++inbr){
// 		dodec0.col(inbr) = half120cell[ nbrs(0,inbr) ].coeffs();
// 	}
// 	dodec0.col(12) = Vector4d(0,0,0,1);
// 	// cout << dodec0.block(0,0,3,1).norm() << endl;
// 	// dumppoints("nbr_icos0.pdb",dodec0.block(0,0,3,13)/dodec0.block(0,0,3,1).norm());
	
// 	return;

// 	Vector3d cen;
// 	double const side0 = min_cube_side(dodec0.block(0,0,3,12),cen);
// 	BOOST_VERIFY( cen.isApprox(Vector3d(0,0,0)) );
// 	double const area0 = side0*side0*side0;
// 	cout << "min cube side: " << side0 << endl;

// 		// cout << dodec0.transpose() << endl;
// 	double mincube = side0;
// 	Matrix4d bestm4;
// 	for(size_t i = 0; i < 1; ++i){
		
// 		Matrix4d randm4; {
// 			randm4 <<   1, 5, 5, 0,
// 			            2, 3, 8, 0,
// 			            3, 4, 7, 0,
// 		                0, 0, 0, 1;
// 			HouseholderQR<Matrix4d> qr(randm4);
// 			randm4 = qr.householderQ();
// 			// cout << randm4 << endl << endl;
// 			// cout << ( randm4 * Vector4d(0,0,0,1) ).transpose() << endl;
// 			BOOST_VERIFY( Vector4d(0,0,0,1).isApprox(randm4*Vector4d(0,0,0,1)) );
// 		}


// 		Matrix<double,4,13>dodec = randm4 * dodec0;
// 		// cout << dodec.transpose() << endl;
// 		// dumppoints("test.pdb",dodec);

// 		double const side = min_cube_side(dodec.block(0,0,3,12),cen);
// 		BOOST_VERIFY( cen.isApprox(Vector3d(0,0,0)) );
// 		if( mincube > side ){
// 			mincube = side;
// 			bestm4 = randm4;
// 			double const area = side*side*side;
// 			cout << "min cube side: " << side << " vol: " << area << " dodec vol frac: " << area / area0 << endl;
// 		}
// 	}


// }

TEST(hecatonicosachoron,neighbor_identity_rots_init){

	vector<Quaterniond> half120cell = make_half120cell();
	Array<uint8_t,60,12> nbrs = make_half120cell_neighbors(half120cell);

	std::vector<Quaterniond> nbrs0;
	for(size_t inbr = 0; inbr < 12; ++inbr)
		nbrs0.push_back( half120cell[ nbrs(0,inbr) ] );

	// size_t ir = 0, ia = 0;
	size_t id = 0;
	BOOST_FOREACH( Quaterniond d, half120cell ){
		Quaterniond const dinv = d.inverse();
		// Quaterniond const dinvd = dinv*d;
		// std::cout << d << "  |  " << dinv << "  |  "  << dinvd << endl;
		// ofstream out(("120cell_"+str(id)+".pdb").c_str());		
		for(size_t inbr = 0; inbr < 12; ++inbr){
			Quaterniond const q = half120cell[ nbrs(id,inbr) ];
			Quaterniond p = to_half_cell( dinv * q );
			// double norm = 100.0/p.coeffs().block(0,0,3,1).norm();
			// std::cout << d << "  |  " << q << "  |  "  << p << endl;
			if( fabs(p.norm()-1.0) > 0.0000001 ) std::exit(-1);
			// if( p.w() < 0.8  ) continue;
			// dumpatom(out,true,++ia,"QUAT","QUT",'A',++ir,norm*p.x(),norm*p.y(),norm*p.z(),1.0,1.0,"C");
			bool in_nbrs0 = false;
			for(size_t jnbr = 0; jnbr < 12; ++jnbr){
				if( p.isApprox(nbrs0[jnbr]) ){
					in_nbrs0 = true;
					// cout << id << " " << inbr << " " << jnbr << endl;
				}
			}
			ASSERT_TRUE(in_nbrs0);
		}
		// out.close();
		++id;
		// break;
	}
}

TEST(hecatonicosachoron,static_inv_check){
	typedef Map<Quaterniond const> QM;
	double const * h120raw = get_h120<double>();
	double const * h120invraw = get_h120inv<double>();	

	Map< Matrix<double,4,60> const > h120( h120raw );
	Map< Matrix<double,4,60> const > h120inv( h120invraw );	
	Map< Array<uint8_t,12,60> const > nbrs( get_h120_nbrs<uint8_t>() );

	for(int i = 0; i < 60; ++i){
		ASSERT_TRUE( Quaterniond(1,0,0,0).isApprox( h120_cellcen<double>(i)*h120_cellceninv<double>(i) ) );
	}
}

#ifdef CXX11
TEST(hecatonicosachoron,neighbor_identity_rots){
	typedef Map<Quaterniond const> QM;
	double const * h120raw = get_h120<double>();

	vector<Quaterniond> h120ref = make_half120cell();
	Array<uint8_t,60,12> nbrsref = make_half120cell_neighbors(h120ref);

	Map< Matrix<double,4,60> const > h120( h120raw );
	Map< Array<uint8_t,12,60> const > nbrs( get_h120_nbrs<uint8_t>() );

	for(int i = 0; i < 60; ++i){
		ASSERT_TRUE( h120ref[i].isApprox( QM( 4*i + h120raw ) ) );
		ASSERT_TRUE( h120ref[i].isApprox( h120_cellcen<double>(i) ) );
		for(int j = 0; j < 12; ++j){
			ASSERT_EQ( nbrsref(i,j), nbrs(j,i) );
		}
	}

	// cout << Quaterniond(h120.col(0)) << endl;
	// cout << Quaterniond(h120.col(1)) << endl;
	// cout << QM( h120raw+0 ) << endl;
	// cout << QM( h120raw+4 ) << endl;
	// cout << nbrs.col(0).transpose().cast<int>() << endl;
	// cout << nbrs.col(1).transpose().cast<int>() << endl;

	std::vector<QM> nbrs0;
	for(size_t inbr = 0; inbr < 12; ++inbr){
		nbrs0.push_back( QM( h120raw + 4*nbrs(inbr,0) ) );
	}
	// BOOST_FOREACH( Quaterniond q , nbrs0 ) cout << q << endl;

	// size_t ir = 0, ia = 0;
	for(int id = 0; id < 60; ++id){
		// QM d( 4*id + h120raw );
		// Quaterniond const dinv = d.inverse();
		QM dinv = h120_cellceninv<double>(id);
		// Quaterniond const dinvd = dinv*d;
		// std::cout << d << "  |  " << dinv << "  |  "  << dinvd << endl;
		// ofstream out(("120cell_"+str(id)+".pdb").c_str());		
		for(size_t inbr = 0; inbr < 12; ++inbr){
			// QM q( h120raw + 4*nbrs(inbr,id) );
			Quaterniond p = to_half_cell( dinv * h120_cellnbr<double>(id,inbr) );
			// double norm = 100.0/p.coeffs().block(0,0,3,1).norm();
			// std::cout << d << "  |  " << q << "  |  "  << p << endl;
			if( fabs(p.norm()-1.0) > 0.0000001 ) std::exit(-1);
			// if( p.w() < 0.8  ) continue;
			// dumpatom(out,true,++ia,"QUAT","QUT",'A',++ir,norm*p.x(),norm*p.y(),norm*p.z(),1.0,1.0,"C");
			bool in_nbrs0 = false;
			for(size_t jnbr = 0; jnbr < 12; ++jnbr){
				if( p.isApprox(nbrs0[jnbr]) ){
					in_nbrs0 = true;
					// cout << id << " " << inbr << " " << jnbr << endl;
				}
			}
			ASSERT_TRUE(in_nbrs0);
		}
		// out.close();
		// break;
	}
}
#endif

TEST(hecatonicosachoron,cell_bounds){
	// cout << h120_cellcen<double>(0) << endl;
	// cout << cellnbr<double>(0,0) << endl;	
	// cout << (h120_cellcen<double>(0).coeffs()+cellnbr<double>(0,0).coeffs()).normalized().transpose() << endl;		
	// cout << std::setprecision(17) << (h120_cellcen<double>(0).coeffs()-cellnbr<double>(0,0).coeffs()).norm() << endl;
}

TEST(hecatonicosachoron,cell_lookup){
	// HecatonicosachoronMap<> map;
	// HecatonicosachoronMap<>::Params p;
	// for( size_t i = 0; i < 60; ++i){
	// 	size_t cell_index = 999;
	// 	map.value_to_params( h120_cellcen<double>(i).matrix(), 0, p, cell_index ); 
	// 	// cout << p << endl;
	// 	ASSERT_EQ( i, cell_index );
	// 	ASSERT_NEAR( p[0], 0.5, 0.0000001 );
	// 	ASSERT_NEAR( p[1], 0.5, 0.0000001 );
	// 	ASSERT_NEAR( p[2], 0.5, 0.0000001 );
	// }

	NEST<3,Matrix3d,HecatonicosachoronMap> nest;
	// ASSERT_EQ( 0, nest.get_index( nest.set_and_get(0,2), 2 ));

    int MAX_RESL = 3;
	#ifdef SCHEME_BENCHMARK
	MAX_RESL += 2;
	#endif
	
	for(int r = 0; r <= MAX_RESL; ++r){
		for(int i = 0; i < (int)nest.size(r); ++i){
			if( nest.set_state(i,r) ){
				size_t ilookup = nest.get_index( nest.value(), r );
				size_t ci0 = i >> (3*r);
				size_t ciL = ilookup >> (3*r);
				if( ci0 == ciL ) ASSERT_EQ( i , ilookup );
			}
		}
	}

}


TEST(hecatonicosachoron,covering){
	// cout << "QuaternionMap Covrad" << endl;
	int NRES = 8;
	int NITER = 10*1000;
	#ifdef SCHEME_BENCHMARK
		NITER = 1000*1000;
	#endif

	std::mt19937 rng((unsigned int)time(0));
	std::normal_distribution<> gauss;
	std::uniform_real_distribution<> uniform;

	NEST<3,Matrix3d,HecatonicosachoronMap> nest;
	for(int r = 0; r <= NRES; ++r){
		double maxdiff=0, avgdiff=0;
		for(int i = 0; i <= NITER; ++i){
			Eigen::Quaterniond q( fabs(gauss(rng)), gauss(rng), gauss(rng), gauss(rng) );
			q.normalize();
			// cout << q << endl;
			Eigen::Quaterniond qcen( nest.set_and_get( nest.get_index(q.matrix(),r) , r ) );
			maxdiff = std::max(maxdiff,q.angularDistance(qcen));
			avgdiff += q.angularDistance(qcen);
		}
		avgdiff /= NITER;
		size_t count = 0; for(size_t i = 0; i < nest.size(r)/60; ++i) if(nest.set_state(i,r)) ++count;
		count *= 60;
		// size_t count = nest.size(r);
		double volfrac = (double)count*(maxdiff*maxdiff*maxdiff)*4.0/3.0*M_PI / 8.0 / M_PI / M_PI;
		double avgfrac = (double)count*(avgdiff*avgdiff*avgdiff)*4.0/3.0*M_PI / 8.0 / M_PI / M_PI;		
		printf("%2i %16lu %10.5f %10.5f %10.5f %10.5f %10.5f\n", 
			r, count, maxdiff*180.0/M_PI, avgdiff*180.0/M_PI, maxdiff/avgdiff, volfrac, avgfrac );
	}

}

// TEST( hecatonicosachoron, visualize ){

// 	std::mt19937 rng((unsigned int)time(0));
// 	std::normal_distribution<> gauss;
// 	std::uniform_real_distribution<> uniform;

// 	// Quaterniond qrand( gauss(rng), gauss(rng), gauss(rng), gauss(rng) );
// 	Quaterniond qrand( 1,0,0,0 );
// 	qrand.normalize();
// 	Vector3d X = qrand*Vector3d(1,0  ,0  );
// 	Vector3d Y = qrand*Vector3d(0,1.2,0  );
// 	Vector3d Z = qrand*Vector3d(0,0  ,1.4);

// 	NEST<3,Matrix3d,HecatonicosachoronMap> nest;
// 	// size_t beg = 0;
// 	// while(!nest.set_state(beg,10)) beg = std::max<size_t>(uniform(rng)*(nest.size(10)-1000),0);
// 	for(size_t r = 0; r <= 4; ++r){
// 		int N = 2 * std::pow(8,r);
// 		int beg = std::max( 0, (int)nest.size(r)/12 - N/2 );
// 		std::ofstream out(("h120_"+boost::lexical_cast<std::string>(r)+".pdb").c_str());
// 		size_t count1 = 0, count2 = 0;
// 		// cout << r << " " << nest.size(r) << " " << (beg>>(4*(10-r))) << endl;
// 		// continue;
// 		// for(size_t i = beg>>(4*(10-r)); i < nest.size(r); ++i){
// 		for(size_t i = beg; i < nest.size(r); ++i){		
// 			++count1;
// 			if( nest.set_state(i,r) ){
// 				++count2;
// 				if( count1 > N) break;
// 				Matrix3d m = nest.value();
// 				// cout << r << " " << i << " " << q.coeffs().transpose() << endl;
// 				 Vector3d ximg = m * X;
// 				 Vector3d yimg = m * Y;
// 				 Vector3d zimg = m * Z;;
// 				 // cout << r << " " << nest.cell_index(i,r) << endl;
// 				 io::dump_pdb_atom(out, nest.cell_index(i,r)<5?"H":"O" ,60*ximg);
// 				 io::dump_pdb_atom(out, nest.cell_index(i,r)<5?"H":"NI",60*yimg);
// 				 io::dump_pdb_atom(out, nest.cell_index(i,r)<5?"H":"N" ,60*zimg);
// 			}
// 		}
// 		out.close();
// 		cout << r << " " << count2 / (double)count1 << endl;
// 	}

// }





}}}}
