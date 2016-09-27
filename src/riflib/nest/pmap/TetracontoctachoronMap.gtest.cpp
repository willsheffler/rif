#include <gtest/gtest.h>

#include "riflib/numeric/util.hpp"
#include "riflib/io/dump_pdb_atom.hpp"

#include "riflib/nest/pmap/TetracontoctachoronMap.hpp"
#include "riflib/nest/NEST.hpp"

#include <Eigen/Dense>
#include <random>

#include <boost/lexical_cast.hpp>

#include <fstream>

namespace scheme { namespace nest { namespace pmap { namespace test {

using namespace Eigen;
using std::cout;
using std::endl;


TEST( TetracontoctachoronMap , cell_lookup ){
	// HecatonicosachoronMap<> map;
	// HecatonicosachoronMap<>::Params p;
	// for( size_t i = 0; i < 60; ++i){
	// 	size_t cell_index = 999;
	// 	map.value_to_params( cellcen<double>(i).matrix(), 0, p, cell_index ); 
	// 	// cout << p << endl;
	// 	ASSERT_EQ( i, cell_index );
	// 	ASSERT_NEAR( p[0], 0.5, 0.0000001 );
	// 	ASSERT_NEAR( p[1], 0.5, 0.0000001 );
	// 	ASSERT_NEAR( p[2], 0.5, 0.0000001 );
	// }

	NEST<3,Matrix3d,TetracontoctachoronMap> nest;
	// ASSERT_EQ( 0, nest.get_index( nest.set_and_get(0,2), 2 ));

	for(int r = 0; r <= 4; ++r){
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

TEST( TetracontoctachoronMap , nside_cell_lookup ){
	int MAX_NSIDE = 11;
	#ifdef SCHEME_BENCHMARK
	MAX_NSIDE = 23;
	#endif
	for(int nside = 1; nside <= MAX_NSIDE; ++nside){
		NEST<3,Matrix3d,TetracontoctachoronMap> nest(nside);
		uint64_t nfail = 0;
		for(int i = 0; i < (int)nest.size(0); ++i){
			if( nest.set_state(i,0) ){
				size_t ilookup = nest.get_index( nest.value(), 0 );
				nfail += ( i != ilookup );
			}
		}
		double frac_fail = (double)nfail / (double)nest.size(0);
		// cout << nside << " " << frac_fail << endl;
		ASSERT_LE( frac_fail, 0.15 );
	}

}


TEST( TetracontoctachoronMap , covering ){
	// cout << "QuaternionMap Covrad" << endl;
	int NRES = 5;
	int NITER = 10*1000;
	#ifdef SCHEME_BENCHMARK
		NITER *= 50;
	#endif

	// std::mt19937 rng((unsigned int)time(0));
	std::mt19937 rng(0);
	std::normal_distribution<> gauss;
	std::uniform_real_distribution<> uniform;

	NEST<3,Matrix3d,TetracontoctachoronMap> nest;

	cout << "resl              N     covrad     avgrad    cov/avg  overcover     avgcov  fracused" << endl;
	for(int r = 0; r <= NRES; ++r){


		double maxdiff=0, avgdiff=0;
		for(int i = 0; i <= NITER; ++i){
			Eigen::Quaterniond q( fabs(gauss(rng)), gauss(rng), gauss(rng), gauss(rng) );
			q.normalize();
			size_t index = nest.get_index(q.matrix(),r);
			// cout << "index: " << index << " sample: " << q.coeffs().transpose() << endl;
			ASSERT_TRUE( nest.set_state( index , r ) );
			Eigen::Quaterniond qcen( nest.value() );
			maxdiff = std::max(maxdiff,q.angularDistance(qcen));
			avgdiff += q.angularDistance(qcen);
		}
		avgdiff /= NITER;
		// size_t count = nest.size(r);
		size_t count = 0; for(size_t i = 0; i < nest.size(r)/24; ++i) if(nest.set_state(i,r)) ++count;
		count *= 24;
		double volfrac = (double)count*(maxdiff*maxdiff*maxdiff)*4.0/3.0*M_PI / 8.0 / M_PI / M_PI;
		double avgfrac = (double)count*(avgdiff*avgdiff*avgdiff)*4.0/3.0*M_PI / 8.0 / M_PI / M_PI;		
		printf("%2i %16lu %10.5f %10.5f %10.5f %10.5f %10.5f %7.5f\n", 
			r,
			count,
			maxdiff*180.0/M_PI,
			avgdiff*180.0/M_PI,
			maxdiff/avgdiff,
			volfrac,
			avgfrac,
			(double)count / nest.size(r)
		);
	}

}

TEST( TetracontoctachoronMap , nside_covering ){
	// cout << "QuaternionMap Covrad" << endl;
	int MAX_NSIDE = 10;
	int NITER = 10*1000;
	#ifdef SCHEME_BENCHMARK
		MAX_NSIDE = 25;
		NITER *= 50;
	#endif

	// std::mt19937 rng((unsigned int)time(0));
	std::mt19937 rng(0);
	std::normal_distribution<> gauss;
	std::uniform_real_distribution<> uniform;


	cout << "nside              N     covrad     avgrad    cov/avg  overcover     avgcov  fracused" << endl;
	for(int nside = 1; nside <= MAX_NSIDE; ++nside){
		NEST<3,Matrix3d,TetracontoctachoronMap> nest(nside);

		double maxdiff=0, avgdiff=0;
		for(int i = 0; i <= NITER; ++i){
			Eigen::Quaterniond q( fabs(gauss(rng)), gauss(rng), gauss(rng), gauss(rng) );
			q.normalize();
			size_t index = nest.get_index(q.matrix(),0);
			// cout << "index: " << index << " sample: " << q.coeffs().transpose() << endl;
			ASSERT_TRUE( nest.set_state( index , 0 ) );
			Eigen::Quaterniond qcen( nest.value() );
			maxdiff = std::max(maxdiff,q.angularDistance(qcen));
			avgdiff += q.angularDistance(qcen);
		}
		avgdiff /= NITER;
		// size_t count = nest.size(0);
		size_t count = 0; for(size_t i = 0; i < nest.size(0)/24; ++i) if(nest.set_state(i,0)) ++count;
		count *= 24;
		double volfrac = (double)count*(maxdiff*maxdiff*maxdiff)*4.0/3.0*M_PI / 8.0 / M_PI / M_PI;
		double avgfrac = (double)count*(avgdiff*avgdiff*avgdiff)*4.0/3.0*M_PI / 8.0 / M_PI / M_PI;		
		printf("%3i %16lu %10.5f %10.5f %10.5f %10.5f %10.5f %7.5f\n", 
			nside,
			count,
			maxdiff*180.0/M_PI,
			avgdiff*180.0/M_PI,
			maxdiff/avgdiff,
			volfrac,
			avgfrac,
			(double)count / nest.size(0)
		);
		std::cout.flush();
	}

}

// TEST( TetracontoctachoronMap , visualize ){

// 	std::mt19937 rng((unsigned int)time(0));
// 	std::normal_distribution<> gauss;
// 	std::uniform_real_distribution<> uniform;

// 	Quaterniond qrand( gauss(rng), gauss(rng), gauss(rng), gauss(rng) );
// 	// Quaterniond qrand( 1,0,0,0 );
// 	qrand.normalize();
// 	Vector3d X = qrand*Vector3d(1,0  ,0  );
// 	Vector3d Y = qrand*Vector3d(0,1.2,0  );
// 	Vector3d Z = qrand*Vector3d(0,0  ,1.4);

// 	NEST<3,Matrix3d,TetracontoctachoronMap> nest;
// 	// size_t beg = 0;
// 	// while(!nest.set_state(beg,10)) beg = std::max<size_t>(uniform(rng)*(nest.size(10)-1000),0);

// 	for(size_t r = 0; r <= 8; ++r){
// 		int N = 8*8*8;
// 		// int beg = std::max( 0, (int)nest.size(r)/12 - N/2 );
// 		int beg = 0;
// 		std::ofstream out(("tcoc_"+boost::lexical_cast<std::string>(r)+".pdb").c_str());
// 		io::dump_pdb_atom(out,  "Z" ,0,Vector3d(0,0,0));
// 		int count1 = 0, count2 = 0;
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
// 				Vector3d ximg = m * X;
// 				Vector3d yimg = m * Y;
// 				Vector3d zimg = m * Z;;
// 				// cout << r << " " << nest.cell_index(i,r) << endl;
// 				// out << "MODEL" << std::endl;
// 				// io::dump_pdb_atom(out,  "Z" ,count2,Vector3d(0,0,0));
// 				io::dump_pdb_atom(out, "O" ,count2,50*ximg);
// 				io::dump_pdb_atom(out, "NI",count2,50*yimg);
// 				io::dump_pdb_atom(out, "N" ,count2,50*zimg);
// 				// out << "ENDMDL" << std::endl;				 
// 			}
// 		}
// 		out.close();
// 		cout << r << " " << count2 / (double)count1 << endl;
// 	}

// }

TEST( TetracontoctachoronMap , check_unit ){

	std::mt19937 rng((unsigned int)time(0));
	std::normal_distribution<> gauss;
	std::uniform_real_distribution<> uniform;

	TetracontoctachoronMap<> map;

	scheme::util::SimpleArray<24,int> cellcount(0);

	// std::ofstream out("bt24_unit.pdb");
	Vector3d lb(9e9,9e9,9e9), ub(-9e9,-9e9,-9e9);
	double lbcorner = 9e9, ubcorner = -9e9;
	for(int i = 0; i < 20000; ++i){

		Quaterniond qrand( gauss(rng), gauss(rng), gauss(rng), gauss(rng) );
		qrand.normalize();

		uint64_t cell_index;
		TetracontoctachoronMap<>::Params params(12345);
		map.value_to_params( qrand.matrix(), 0, params, cell_index );
		cellcount[cell_index]++;
		if( cell_index != 0 ) continue;

		Vector3d p(params[0],params[1],params[2]);
		lb[0] = fmin(lb[0],p[0]);
		lb[1] = fmin(lb[1],p[1]);
		lb[2] = fmin(lb[2],p[2]);
		ub[0] = fmax(ub[0],p[0]);
		ub[1] = fmax(ub[1],p[1]);
		ub[2] = fmax(ub[2],p[2]);
		lbcorner = fmin( lbcorner, p[0]+p[1]+p[2] );
		ubcorner = fmax( ubcorner, p[0]+p[1]+p[2] );		

		// io::dump_pdb_atom(out, "O" ,i/1000, 100*p);

	}
	// out.close();
	// cout << cellcount << endl;

	// cout << "bt24 unit lb: " << lb.transpose() << endl;
	// cout << "bt24 unit ub: " << ub.transpose() << endl;	
	ASSERT_LE( lb[0], 0.03 );
	ASSERT_LE( lb[1], 0.03 );
	ASSERT_LE( lb[2], 0.03 );		
	ASSERT_GE( lb[0], 0.00 );
	ASSERT_GE( lb[1], 0.00 );
	ASSERT_GE( lb[2], 0.00 );		
	ASSERT_LE( ub[0], 1.00 );
	ASSERT_LE( ub[1], 1.00 );
	ASSERT_LE( ub[2], 1.00 );		
	ASSERT_GE( ub[0], 0.97 );
	ASSERT_GE( ub[1], 0.97 );
	ASSERT_GE( ub[2], 0.97 );

	EXPECT_GE( ubcorner , 1.5 + 0.5/(sqrt(2.0)-1.0) -0.2 );
	EXPECT_LE( ubcorner , 1.5 + 0.5/(sqrt(2.0)-1.0) );
	EXPECT_LE( lbcorner , 1.5 - 0.5/(sqrt(2.0)-1.0) +0.2 );
	EXPECT_GE( lbcorner , 1.5 - 0.5/(sqrt(2.0)-1.0) );
	// ASSERT_GE( ubcorner ,  );

}

// nside              N     covrad     avgrad    cov/avg  overcover     avgcov  fracused
//   1               24   62.76235   40.72950    1.54096    1.67355    0.45737 1.00000
//   2              192   39.38715   20.63925    1.90836    3.30899    0.47612 1.00000
//   3              648   26.80538   13.82587    1.93878    3.52024    0.48304 1.00000
//   4             1536   20.18428   10.38848    1.94295    3.56256    0.48571 1.00000
//   5             3000   16.18733    8.31389    1.94702    3.58904    0.48626 1.00000
//   6             5184   13.46439    6.92945    1.94307    3.56908    0.48651 1.00000
//   7             8232   11.58217    5.94196    1.94922    3.60750    0.48711 1.00000
//   8            12288   10.09158    5.19889    1.94111    3.56197    0.48702 1.00000
//   9            17496    9.00114    4.62175    1.94756    3.59884    0.48718 1.00000
//  10            24000    8.08815    4.16021    1.94417    3.58170    0.48740 1.00000

}}}}
