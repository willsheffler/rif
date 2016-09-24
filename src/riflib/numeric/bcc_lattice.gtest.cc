#include <gtest/gtest.h>

#include "riflib/numeric/bcc_lattice.hh"
#include "riflib/io/dump_pdb_atom.hh"
#include "riflib/util/Timer.hh"
#include <fstream>
#include <random>
#include <boost/foreach.hpp>
#include <iterator>     // std::back_inserter
#include <boost/format.hpp>
#include <sparsehash/dense_hash_set>

#include <Eigen/Geometry>

namespace scheme { namespace numeric { namespace test {

using std::cout;
using std::endl;

// TEST(TEMPORARY,gm_20140905){
// 	int NSAMP = 1000000;
// 	std::mt19937 r((unsigned int)time(0));
// 	std::uniform_real_distribution<> u;
// 	Eigen::Matrix<double,6,6> dis; dis.fill(0);
// 	for(int idata = 1; idata <= 6; ++idata ){
// 		for(int irep  = 1; irep <= idata; ++irep ){
// 			int totbin = 16<<(idata-1)*2;
// 			int nbins = std::pow( totbin, 1.0/irep );
// 			cout << idata << " " << irep << " " << totbin << " " << std::pow(nbins,irep) << endl;
// 			for(int iter = 0; iter < NSAMP; ++iter){
// 				util::SimpleArray<6,double> samp(0.5);
// 				for(int i=0; i<idata; ++i) samp[i] = u(r);
// 				util::SimpleArray<6,double> rep(0.5);
// 				for(int i=0; i<irep; ++i){
// 					rep[i] = ((int)(samp[i]*nbins)+u(r))/nbins;
// 				}
// 				// cout << idata << " " << irep << " " << samp << endl;
// 				// cout << idata << " " << irep << " " << rep << endl;
// 				// cout << endl;
// 				// std::exit(-1);
// 				dis(irep-1,idata-1) += (samp-rep).norm();
// 			}
// 		}
// 	}
// 	dis = dis / (double)NSAMP * 100;
// 	cout << dis << endl;
// }

TEST(bcc_lattice,centers_map){
	typedef util::SimpleArray<3,double> V;
	typedef util::SimpleArray<3,uint64_t> I;

	// {
	// 	BCC<3,double> bcc(I(2,2,2),V(0,0,0),V(1,1,1));
	// 	for(int i = 0; i < bcc.size(); ++i)
	// 		cout << i << "\t" << bcc[i] << endl;

	// }

	BCC<3,double> bcc(I(3,5,7),V(0,0,0),V(6,10,14));

	// std::ofstream out("test.pdb");
	for(size_t i = 0; i < bcc.size(); i+=2){
		ASSERT_EQ( i, bcc[bcc[i]] );
		// printf("%6lu %10.6f %10.6f %10.6f\n",i,bcc[i][0],bcc[i][1],bcc[i][2]);
		// io::dump_pdb_atom(out,i,bcc[i]*3.0);
		// cout << i << " " << bcc[i] << endl;
	}
	// out.close();

}

	// double const R3approx = 0.558099;
	// double const R4approx = 0.701687;
	// double const R5approx = 0.742306;
	// double const R6approx = 0.845359;
	// double const R7approx = 0.882879;
template<int N, class F, class S>
F
test_bcc_performance(
	size_t NSAMP,
	S const Nside,
	F const Width
){
	typedef util::SimpleArray<N,F> V;
	typedef util::SimpleArray<N,S> I;
	std::mt19937 rng((unsigned int)time(0));
	std::uniform_real_distribution<> runif;
	BCC<N,F,S> bcc(I(Nside),V(-Width/2),V(Width/2));

	std::vector<V> samples(NSAMP);
	for(int i = 0; i < NSAMP; ++i)
		for(int j = 0; j < N; ++j)
			// samples[i][j] = runif(rng)*9.4+0.3;
			samples[i][j] = runif(rng)*9.0-4.5;			

	std::vector<size_t> indices(NSAMP);
	util::Timer<> lookup_time;
	for(int i = 0; i < NSAMP; ++i)
		indices[i] = bcc[samples[i]];
	cout << "BCC DIM " << N << " lookup rate: " << (double)NSAMP / lookup_time.elapsed() << " sec / ";

	std::vector<V> centers(NSAMP);
	util::Timer<> getval_time;
	for(int i = 0; i < NSAMP; ++i)
		centers[i] = bcc[indices[i]];
	cout << " getval rate: " << (double)NSAMP / getval_time.elapsed() << " sec" << endl;


	F maxdiff = 0;
	for(int i = 0; i < NSAMP; ++i)
		maxdiff = std::max( maxdiff, (samples[i]-centers[i]).squaredNorm() );

	maxdiff = sqrt(maxdiff);
	F frac = bcc.width_[0] * sqrt(N)/2.0 / maxdiff;
	F improvement = 1.0; for(int i = 0; i < N; ++i) improvement *= frac;
	cout << "     improvement over cubic: " << improvement / 2.0 << " cov: " << maxdiff / Width * Nside
	     << " vs. " << sqrt(N)/2.0 << endl; // 2 x num samp as cubic

	return maxdiff;
}

TEST(bcc_lattice,performance){
	size_t NITER = 50*1000;
	#ifdef SCHEME_BENCHMARK
	NITER *= 50;
	#endif
	size_t Nside = 100;
	double Width = 10.0;
	double const R3test = test_bcc_performance<3,double,uint64_t>( NITER, Nside, Width );
	                      test_bcc_performance<4,double,uint64_t>( NITER, Nside, Width );
	                      test_bcc_performance<5,double,uint64_t>( NITER, Nside, Width );
	                      test_bcc_performance<6,double,uint64_t>( NITER, Nside, Width );
	                      test_bcc_performance<7,double,uint64_t>( NITER, Nside, Width );
	double const R3 = std::pow(2.0,-5.0/3.0)*sqrt(5) / std::pow(2.0,1.0/3.0) * Width / Nside;
	ASSERT_LE( R3test     , R3 );
	ASSERT_GT( R3test*1.1 , R3 );	
}


template<int N, class F, class S>
F
test_bcc_inradius(){
	typedef util::SimpleArray<N,F> V;
	typedef util::SimpleArray<N,S> I;
	std::mt19937 rng((unsigned int)time(0));
	std::uniform_real_distribution<> runif;
	std::normal_distribution<> rnorm;	
	S const Nside = 5;
	BCC<N,F,S> bcc(I(5),V(-(F)Nside/2.0),V((F)Nside/2.0));
	BOOST_VERIFY( bcc[bcc[V(0.0)]] == V(0.0) );
	BOOST_VERIFY( bcc.width_ == V(1.0) );
	S const i0 = bcc[V(0)];

	double const RNapprox[5] = { 0.558099, 0.701687, 0.742306, 0.845359, 0.882879 };
	double const Rapprox = RNapprox[N-3];

	int NSAMP = 50*1000;

	double min_inrad = Rapprox;
	for(int i = 0; i < NSAMP; ++i){
		V samp;
		for(int j = 0; j < N; ++j) samp[j] = rnorm(rng);
		samp.normalize();
		double const radius = runif(rng)*min_inrad*0.1 + 0.9*min_inrad;
		samp *= radius;
		if( bcc[samp] != i0 ){
			min_inrad = radius;
		}

	}
	return min_inrad;
}

TEST(bcc_lattice,inradius){
	ASSERT_NEAR( (test_bcc_inradius<3,double,uint64_t>()), 0.433015, 0.03 );
	ASSERT_NEAR( (test_bcc_inradius<4,double,uint64_t>()), 0.500000, 0.03 );
	ASSERT_NEAR( (test_bcc_inradius<5,double,uint64_t>()), 0.500000, 0.03 );
	ASSERT_NEAR( (test_bcc_inradius<6,double,uint64_t>()), 0.500000, 0.03 );
	ASSERT_NEAR( (test_bcc_inradius<7,double,uint64_t>()), 0.500000, 0.03 );
}

template<int N, class F, class S> 
F 
test_bcc_neighbors( size_t NSAMP ){
	typedef util::SimpleArray<N,F> V;
	typedef util::SimpleArray<N,S> I;
	std::mt19937 rng((unsigned int)time(0));
	std::uniform_real_distribution<> runif;
	std::normal_distribution<> rnorm;	
	S const Nside = 5;
	BCC<N,F,S> bcc(I(5),V(-(F)Nside/2.0),V((F)Nside/2.0));
	BOOST_VERIFY( bcc[bcc[V(0.0)]] == V(0.0) );
	BOOST_VERIFY( bcc.width_ == V(1.0) );
	S const i0 = bcc[V(0)];
	std::vector<size_t> nbrs,nbrs_we;		
	bcc.neighbors( i0, std::back_inserter(nbrs   ), false );
	bcc.neighbors( i0, std::back_inserter(nbrs_we), true );		

	Cubic<N,F,S> cubic(I(5),V(-(F)Nside/2.0),V((F)Nside/2.0));
	BOOST_VERIFY( cubic[cubic[V(0.0)]] == V(0.0) );
	BOOST_VERIFY( cubic.width_ == V(1.0) );
	S const i0cubic = cubic[V(0)];
	std::vector<size_t> nbrs_cubic;		
	cubic.neighbors( i0cubic, std::back_inserter(nbrs_cubic) );

	/////////////////////////////////////////////////////
	// test neighbor coverage
	///////////////////////////////////////////////////////////////

	F maxrad_99 = 0;
	// F const RNapprox[5] = { 0.558099, 0.701687, 0.742306, 0.845359, 0.882879 };
	// F const Rapprox = RNapprox[N-3];
	F inrad = 0.5;
	if(N==3) inrad = 0.433015;
	// F const radius = 1.9*inrad;
	for(F radius = 0; radius < 3.0*inrad; radius += inrad/10.0){
	// F radius = 2.1*inrad; {

		int sum_in_nbrs=0, sum_in_nbrs_we=0;
		for(int i = 0; i < NSAMP; ++i){

			// pick random point in i0 cell
			V cen; for(int j = 0; j < N; ++j) cen[j] = (runif(rng)-0.5)*2.0*inrad;
			if( bcc[cen] != i0 ){ --i; continue; }

			V samp;
			for(int j = 0; j < N; ++j) samp[j] = rnorm(rng);
			samp *= (radius/samp.norm());
			S index = bcc[samp+cen];
			// if( index == i0 ){ --i; continue; }
			if( std::find( nbrs_we.begin(), nbrs_we.end(), index ) != nbrs_we.end() ){
				++sum_in_nbrs_we;
				if( std::find( nbrs.begin(), nbrs.end(), index ) != nbrs.end() ){
					++sum_in_nbrs;
				}
			}
		}
		if( (F)sum_in_nbrs_we/NSAMP > 0.99 ) maxrad_99 = radius;

		// uncomment for print report
		if( fabs( radius - 2.0*inrad ) > 0.001 ) continue;
		int sum_in_nbrs_cubic=0;
		for(int i = 0; i < NSAMP; ++i){
			// pick random point in i0 cell
			V cen; for(int j = 0; j < N; ++j) cen[j] = (runif(rng)-0.5);
			BOOST_VERIFY( cubic[cen] == i0cubic );

			V samp;
			for(int j = 0; j < N; ++j) samp[j] = rnorm(rng);
			samp *= (radius/inrad*0.5/samp.norm());
			S index = cubic[samp+cen];
			// if( index == i0 ){ --i; continue; }
			if( std::find( nbrs_cubic.begin(), nbrs_cubic.end(), index ) != nbrs_cubic.end() ){
				++sum_in_nbrs_cubic;
			}
		}
		printf("%i %7.3f %9.7f %9.7f %9.7f\n",
			N,
			radius/inrad,
			(F)sum_in_nbrs/NSAMP,
			(F)sum_in_nbrs_we/NSAMP,
			(F)sum_in_nbrs_cubic/NSAMP			
		);

	}

	// ///////////////////////////////////////////////////////////////////
	// // // dump neighbors
	// //////////////////////////////////////////////////////////////////
	// if(N != 3) return;
	// #define LAT bcc
	// std::ofstream out_bcc("bcc.pdb");
	// for(int i = 0; i < LAT.size(); ++i) io::dump_pdb_atom(out_bcc,i,LAT[i]*10.0);
	// out_bcc.close();
	// for(int i = 0; i < LAT.size(); ++i){
	// 	std::vector<size_t> nbrs;
	// 	LAT.neighbors( i, std::back_inserter(nbrs), true );
	// 	// BOOST_FOREACH(size_t nbr, nbrs) cout << nbr << " " << LAT[nbr] << endl;
	// 	std::string s = boost::str(boost::format("%4i") % i);
	// 	std::ofstream out("test_"+s+".pdb");
	// 	BOOST_FOREACH(size_t nbr, nbrs) io::dump_pdb_atom(out,nbr,LAT[nbr]*10.0);
	// 	out.close();
	// }

	return maxrad_99;
}

TEST(bcc_lattice,neighbors){
	size_t NITER = 1000;
	#ifdef SCHEME_BENCHMARK
	NITER *= 50;
	#endif
	// for(int i = 3; i < 8; ++i){
	// 	int nc = std::pow(3,i);
	// 	int nbccFC = 1+2*i+std::pow(2,i);
	// 	int nbccFCE = nbccFC + i*(i-1)/2 * 4;
	// 	cout << "Nnbrs: " << i <<" cubic "<< nc << " bccFC " << nbccFC << " bccFCE " << nbccFCE << endl;
	// }
	// Nnbrs: 3 cubic   27 bccFC  15 bccFCE  27
	// Nnbrs: 4 cubic   81 bccFC  25 bccFCE  49
	// Nnbrs: 5 cubic  243 bccFC  43 bccFCE  83
	// Nnbrs: 6 cubic  729 bccFC  77 bccFCE 137
	// Nnbrs: 7 cubic 2187 bccFC 143 bccFCE 227
	double v1 = test_bcc_neighbors<3,double,uint64_t>(NITER);
	double v2 = test_bcc_neighbors<4,double,uint64_t>(NITER);
	double v3 = test_bcc_neighbors<5,double,uint64_t>(NITER);
	double v4 = test_bcc_neighbors<6,double,uint64_t>(NITER);
	double v5 = test_bcc_neighbors<7,double,uint64_t>(NITER);
	EXPECT_LE( 0.71, v1 ); // these are total approximations
	EXPECT_LE( 0.67, v2 ); // these are total approximations
	EXPECT_LE( 0.62, v3 ); // these are total approximations
	EXPECT_LE( 0.53, v4 ); // these are total approximations
	EXPECT_LE( 0.49, v5 ); // these are total approximations this one in particular
	// cout << test_bcc_neighbors<3,double,size_t>(1000000) << endl; // 0.779427
	// cout << test_bcc_neighbors<4,double,size_t>(1000000) << endl; // 0.75
	// cout << test_bcc_neighbors<5,double,size_t>(1000000) << endl; // 0.7
	// cout << test_bcc_neighbors<6,double,size_t>(1000000) << endl; // 0.65
	// cout << test_bcc_neighbors<7,double,size_t>(1000000) << endl; // 0.6
}

// TEST(bcc_lattice,coverage_transform_7d){
// 	using namespace Eigen;
// 	typedef Transform<double,3,AffineCompact> Xform;
// 	typedef util::SimpleArray<7,double> V;
// 	typedef util::SimpleArray<7,size_t> I;
// 	typedef Matrix<double,7,1> Vector7d;
// 	std::mt19937 mt((unsigned int)time(0));
// 	std::normal_distribution<> rnorm;
// 	size_t Nside = 64;
// 	V bounds = V(1.5,1.5,1.5,1.5,10,10,10);
// 	BCC<7,double> bcc(I(Nside),-bounds,bounds);

// 	Matrix<double,3,6> pts0;
// 	pts0 <<  1, 0, 0,-2, 0, 0,
// 	         0, 1, 0, 0,-1, 0,
// 	         0, 0, 1, 0, 0,-1;
// 	pts0.colwise() -= pts0.rowwise().sum()/pts0.cols();

// 	Xform X( AngleAxisd(2,Vector3d(1,2,3)) );

// 	// cout << pts0 << endl;
// 	// cout << X*pts0 << endl;

// 	int const NSAMP = 3;
// 	for(int i = 0; i < NSAMP; ++i){
// 		Vector4d quat( rnorm(mt), rnorm(mt), rnorm(mt), rnorm(mt) ); quat.normalize();
// 		Vector3d trans( rnorm(mt), rnorm(mt), rnorm(mt) );
// 		Vector7d samp;
// 		samp.block(0,0,4,1) = quat;
// 		samp.block(4,0,3,1) = trans;
// 		// cout << samp.transpose() << endl;
// 		V tmp = bcc[ bcc[ samp ] ];
// 		Vector7d cen;
// 		for(int i = 0; i < 7; ++i) cen[i] = tmp[i];

// 		cout << samp.transpose() << endl;
// 		cout << cen.transpose() << endl;
// 		cout << endl;

// // 	std::vector<size_t> indices(NSAMP);
// // 	util::Timer<> lookup_time;
// // 	for(int i = 0; i < NSAMP; ++i)
// // 		indices[i] = bcc[samples[i]];
// // 	cout << N << " lookup: " << (double)NSAMP / lookup_time.elapsed() << " sec" << endl;

// // 	std::vector<V> centers(NSAMP);
// // 	util::Timer<> getval_time;
// // 	for(int i = 0; i < NSAMP; ++i)
// // 		centers[i] = bcc[indices[i]];

// 	}


// }



template<int N, class F, class S> 
F 
test_bcc_children( size_t NSAMP ){
	typedef util::SimpleArray<N,F> V;
	typedef util::SimpleArray<N,S> I;
	std::mt19937 rng((unsigned int)time(0));
	std::uniform_real_distribution<> runif;
	std::normal_distribution<> rnorm;	
	S const Nside = 5;
	BCC<N,F,S> bcc_parent(I(5),V(-(F)Nside),V((F)Nside));
	BCC<N,F,S> bcc(I(5),V(-(F)Nside/2.0),V((F)Nside/2.0));
	BOOST_VERIFY( bcc[bcc[V(0.0)]] == V(0.0) );
	BOOST_VERIFY( bcc.width_ == V(1.0) );
	S const i0 = bcc[V(0)];
	std::vector<size_t> nbrs,nbrs_we;		
	bcc.neighbors( i0, std::back_inserter(nbrs   ), false );
	bcc.neighbors( i0, std::back_inserter(nbrs_we), true );		

	/////////////////////////////////////////////////////
	// test neighbor coverage
	///////////////////////////////////////////////////////////////

	// F const RNapprox[5] = { 0.558099, 0.701687, 0.742306, 0.845359, 0.882879 };
	// F const Rapprox = RNapprox[N-3];
	F inrad = 0.5;
	if(N==3) inrad = 0.433015;

	int sum_in_nbrs=0, sum_in_nbrs_we=0;
	for(int i = 0; i < NSAMP; ++i){

		// pick random point in i0 parent cell
		V samp; for(int j = 0; j < N; ++j)
			samp[j] = (runif(rng)-0.5)*4.0*inrad;
		if( bcc_parent[samp] != i0 ){ --i; continue; }

		S index = bcc[samp];
		// if( index == i0 ){ --i; continue; }
		if( std::find( nbrs_we.begin(), nbrs_we.end(), index ) != nbrs_we.end() ){
			++sum_in_nbrs_we;
			if( std::find( nbrs.begin(), nbrs.end(), index ) != nbrs.end() ){
				++sum_in_nbrs;
			}
		}
	}

	//int nc = std::pow(3,N);
	int nbccFC = 1+2*N+std::pow(2,N);
	int nbccFCE = nbccFC + N*(N-1)/2 * 4;
	printf("BCC child coverage: %i %3i %9.7f   %3i %9.7f\n",
		N, 
		nbccFC, (F)sum_in_nbrs/NSAMP,
		nbccFCE, (F)sum_in_nbrs_we/NSAMP
	);
	return 0;
}

TEST(bcc_lattice,children){
	int NITER = 5*1000;
	#ifdef SCHEME_BENCHMARK
	NITER *= 30;
	#endif

	cout << "BCC               DIM Nfc   frac_fc  Nfce  frac_fce" << std::endl;
	test_bcc_children<3,double,uint64_t>( NITER );
	test_bcc_children<4,double,uint64_t>( NITER );
	test_bcc_children<5,double,uint64_t>( NITER );
	test_bcc_children<6,double,uint64_t>( NITER );
	test_bcc_children<7,double,uint64_t>( NITER );
}


}}}

