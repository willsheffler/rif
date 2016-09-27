#include <gtest/gtest.h>

#include "riflib/numeric/util.hpp"
#include "riflib/io/dump_pdb_atom.hpp"

#include "riflib/numeric/geom_4d.hpp"

#include <Eigen/Dense>
#include <random>

#include "riflib/util/Timer.hpp"

#include <fstream>

namespace scheme { namespace nest { namespace pmap {

using namespace Eigen;
using std::cout;
using std::endl;

TEST( geom_4d, quat_half_cell ){
	ASSERT_EQ( -0.0, 0.0 );
	ASSERT_EQ( numeric::to_half_cell(Eigen::Quaterniond( -1, 1, 1, 1 ) ).w(), 1.0 );
	ASSERT_EQ( numeric::to_half_cell(Eigen::Quaterniond(  0,-1, 1, 1 ) ).x(), 1.0 );
	ASSERT_EQ( numeric::to_half_cell(Eigen::Quaterniond(  0, 0,-1, 1 ) ).y(), 1.0 );
	ASSERT_EQ( numeric::to_half_cell(Eigen::Quaterniond(  0, 0, 0,-1 ) ).z(), 1.0 );
}



TEST( geom_4d , tetracontoctachoron_cell_lookup )
{
	typedef double Float;
	typedef uint64_t Index;
	typedef Eigen::Matrix<Float,4,1> V4;

	Map<Matrix<Float,48,4,RowMajor>const> t24( numeric::get_raw_48cell<Float>() );
	// ASSERT_EQ( V4(0,0,0,-1), t24.row(7).transpose() );
	// for(int i = 0; i < 48; ++i) ASSERT_FLOAT_EQ( t24.row(i).norm(), 1.0 );

	// std::ofstream out("test.pdb");
	// for(int i = 0; i < 48; ++i){
	// 	if( t24.row(i).block(0,0,1,3).norm() < 0.9 )
	// 		io::dump_pdb_atom(out,i,10*t24.row(i));
	// }
	// out.close();

	std::mt19937 mt((unsigned int)time(0));
	std::normal_distribution<> rnorm;
	std::uniform_real_distribution<> runif;

	int NITER = 200*1000;
	#ifdef SCHEME_BENCHMARK
		NITER *= 50;
	#endif

	std::vector<V4> samp(NITER);
	std::vector<Index> cell(NITER),cell2(NITER);

	for(int i = 0; i < NITER; ++i){
		V4 quat(rnorm(mt),rnorm(mt),rnorm(mt),rnorm(mt));
		samp[i] = quat.normalized();
	}
		// V4 const quat_pos = quat.cwiseAbs();

	util::Timer<> naive;
	for(int i = 0; i < NITER; ++i){
		// (t24*samp[i]).cwiseAbs().maxCoeff(&cell[i]);
		(t24*samp[i]).maxCoeff(&cell[i]);		
	}
	cout << "bt24 naive rate:  " << (Float)NITER / naive.elapsed_nano() << endl;

	util::Timer<> clever;
	for(int i = 0; i < NITER; ++i){
		numeric::get_cell_48cell( samp[i], cell2[i] );
		// this is slower !?!
		// Float mx = std::max(std::max(hyperface_dist,corner_dist),edge_dist);
		// cell2[i] = hyperface_dist==mx ? facecell : (corner_dist==mx ? cornercell+8 : edgecell+24);


	}
	cout << "bt24 clever rate: " << (double)NITER / clever.elapsed_nano() << endl;

	for(int i = 0; i < NITER; ++i){
		if( cell[i] != cell2[i] ){
			ASSERT_FLOAT_EQ(
			     t24.row(cell [i]).dot(samp[i]) ,
			     t24.row(cell2[i]).dot(samp[i]) );
		} else { 
			ASSERT_EQ( cell[i], cell2[i] );
		}
	}

}

TEST( geom_4d , tetracontoctachoron_half_cell_lookup )
{
	typedef double Float;
	typedef uint64_t Index;
	typedef Eigen::Matrix<Float,4,1> V4;

	Map<Matrix<Float,24,4,RowMajor>const> t24h( numeric::get_raw_48cell_half<Float>() );

	for(int i = 0; i < 24; ++i){
		size_t cell;
		V4 tmp = t24h.row(i);
		numeric::get_cell_48cell_half( tmp, cell );
		ASSERT_EQ( cell, i );
		// if( i > 11 ) continue;
		tmp = -tmp;
		numeric::get_cell_48cell_half( tmp, cell );
		ASSERT_EQ( cell, i );
	}

	// ASSERT_EQ( V4(0,0,0,-1), t24h.row(7).transpose() );
	// for(int i = 0; i < 48; ++i) ASSERT_FLOAT_EQ( t24h.row(i).norm(), 1.0 );

	// std::ofstream out("test.pdb");
	// for(int i = 0; i < 24; ++i){
	// 	if( t24h.row(i).block(0,0,1,3).norm() < 0.9 )
	// 		io::dump_pdb_atom(out,i,10*t24h.row(i));
	// }
	// out.close();

	std::mt19937 mt((unsigned int)time(0));
	std::normal_distribution<> rnorm;
	std::uniform_real_distribution<> runif;

	int NITER = 200*1000;
	#ifdef SCHEME_BENCHMARK
		NITER *= 50;
	#endif

	std::vector<V4> samp(NITER);
	std::vector<Index> cell(NITER),cell2(NITER);

	for(int i = 0; i < NITER; ++i){
		V4 quat(rnorm(mt),rnorm(mt),rnorm(mt),rnorm(mt));
		samp[i] = quat.normalized();
	}
		// V4 const quat_pos = quat.cwiseAbs();

	util::Timer<> naive;
	for(int i = 0; i < NITER; ++i){
		// (t24h*samp[i]).cwiseAbs().maxCoeff(&cell[i]);
		(t24h*samp[i]).cwiseAbs().maxCoeff(&cell[i]);		
	}
	cout << "hbt24 naive rate:  " << (Float)NITER / naive.elapsed_nano() << endl;

	util::Timer<> clever;
	for(int i = 0; i < NITER; ++i){
		numeric::get_cell_48cell_half( samp[i], cell2[i] );
	}
	cout << "hbt24 clever rate: " << (double)NITER / clever.elapsed_nano() << endl;

	for(int i = 0; i < NITER; ++i){
		// if( cell[i] > 11 ) continue;
		// if( cell[i] < 12 ) ASSERT_EQ( cell[i], cell2[i] );
		if( cell[i] != cell2[i] ){
			ASSERT_FLOAT_EQ(
			     t24h.row(cell [i]).dot(samp[i]) ,
			     t24h.row(cell2[i]).dot(samp[i]) );
		} else { 
			ASSERT_EQ( cell[i], cell2[i] );
		}
	}

}

}}}
